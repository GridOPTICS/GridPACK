# -------------------------------------------------------------
# file: gridpack/mcp/study.py
# -------------------------------------------------------------
# Studies: one derived case and one job per point, ranked at the end.
# -------------------------------------------------------------

"""Load sweeps and dynamic N-1 screens over the job runner."""

import time
import uuid

from .events import raw_elements
from . import post

__all__ = ["KINDS", "Study", "n1_branches", "rank"]

KINDS = ("load_sweep", "dynamic_n1")


class Study:
    """A study's items; each is a derived case and the job run on it."""

    def __init__(self, kind, base, params):
        self.id = "s" + uuid.uuid4().hex[:8]
        self.kind = kind
        self.base = base
        self.params = params
        self.items = []            # {"name", "case", "job_id", ...point}
        self.created = time.time()

    def status(self, runner):
        states = {}
        jobs = []
        for it in self.items:
            job = runner.get(it["job_id"])
            states[job.state] = states.get(job.state, 0) + 1
            jobs.append({"name": it["name"], "case": it["case"],
                         "job_id": it["job_id"], "state": job.state})
        return {"id": self.id, "kind": self.kind, "base": self.base,
                "params": self.params, "items": len(self.items),
                "states": states,
                "pending": states.get("queued", 0) + states.get("running", 0),
                "jobs": jobs}


def n1_branches(raw, only=None):
    """``[(from, to, ckt)]`` to outage: every branch and transformer in the
    RAW, or those named in ``only`` as "from-to" or "from-to ckt"."""
    _, branches, _ = raw_elements(raw)
    allb = sorted(branches)
    if not only:
        return allb
    out = []
    for spec in only:
        parts = str(spec).replace("-", " ").split()
        try:
            i, j = int(parts[0]), int(parts[1])
        except (IndexError, ValueError):
            raise ValueError("branch %r is not \"from-to\" or \"from-to ckt\""
                             % spec)
        ckt = parts[2] if len(parts) > 2 else None
        hits = [b for b in allb if {b[0], b[1]} == {i, j}
                and (ckt is None or b[2] == ckt)]
        if not hits:
            raise ValueError("no branch %s in the RAW" % spec)
        if len(hits) > 1:
            raise ValueError("branch %d-%d has circuits %s; name one"
                             % (i, j, ", ".join(b[2] for b in hits)))
        out += hits
    return out


def rank(study, runner):
    """The finished items as a table, most severe first, plus failures."""
    rows, failed, pending = [], [], 0
    for it in study.items:
        job = runner.get(it["job_id"])
        if job.state in ("queued", "running"):
            pending += 1
            continue
        result = job.result
        if job.state != "done" or result is None:
            failed.append({"name": it["name"], "job_id": it["job_id"],
                           "error": job.error})
            continue
        row = _sweep_row(it, result) if study.kind == "load_sweep" \
            else _n1_row(it, result, study.params)
        row.update(name=it["name"], case=it["case"], job_id=it["job_id"])
        rows.append(row)
    if study.kind == "load_sweep":
        rows.sort(key=lambda r: r["factor"])
    else:
        rows.sort(key=lambda r: (not r["unstable"], -r["angle_spread_peak_deg"],
                                 -r["peak_speed_deviation"]))
    return {"id": study.id, "kind": study.kind, "base": study.base,
            "finished": len(rows), "pending": pending, "table": rows,
            "failed": failed}


def _sweep_row(it, r):
    buses, branches = r.get("buses", []), r.get("branches", [])
    row = {"factor": it["factor"], "converged": r.get("converged"),
           "iterations": r.get("iterations")}
    if buses:
        lo = min(buses, key=lambda b: b["voltage"])
        slack = [b for b in buses if b.get("type") == 3]
        row.update(load_mw=round(sum(b["pLoad"] for b in buses), 3),
                   generation_mw=round(sum(b["pGen"] for b in buses), 3),
                   losses_mw=round(sum(b.get("pLoss", 0.0) for b in branches), 3),
                   min_voltage=lo["voltage"], min_voltage_bus=lo["busId"],
                   slack_mw=slack[0]["pGen"] if slack else None,
                   voltage_violations=sum(1 for b in buses
                                          if not 0.9 <= b["voltage"] <= 1.1))
        rated = [b for b in branches if b.get("rateA", 0) > 0]
        row["max_loading_percent"] = max((b["loadingPercent"] for b in rated),
                                         default=None)
    return row


def _n1_row(it, r, params):
    st = post.dynamic_stats(r, "dsf", params.get("fault_start", 0.0)) or {}
    return {"branch": it["branch"],
            "unstable": bool(st.get("loss_of_synchronism") or st.get("runaway_speed")),
            "angle_spread_peak_deg": st.get("angle_spread_peak", 0.0),
            "angle_spread_final_deg": st.get("angle_spread_final", 0.0),
            "peak_speed_deviation": st.get("peak_speed_deviation", 0.0),
            "final_speed_deviation": st.get("final_speed_deviation", 0.0),
            "damping_ratio": st.get("damping_ratio")}
