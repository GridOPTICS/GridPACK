# -------------------------------------------------------------
# file: gridpack/mcp/post.py
# -------------------------------------------------------------
# Post-processing of job results: series, reports, comparisons.
# -------------------------------------------------------------

"""Read a finished job's ``--json`` result back as something a reader can
use: a decimated time series, a short text verdict, or a diff of two runs.
"""

import fnmatch
import math
import xml.etree.ElementTree as ET
from pathlib import Path

__all__ = ["compare", "dynamic_stats", "events_from_xml", "report", "series",
           "table"]

_TIME = {"dsf": "t", "hadrec": "time"}
_TABLE = {"dsf": "watch", "hadrec": "observations"}
_DEG = 180.0 / math.pi
# Voltage band and loading limit of PFResult.violations().
_VMIN, _VMAX, _LOAD = 0.9, 1.1, 100.0


def table(result, app):
    """``(columns, rows)`` of a dynamic run, or None for other apps."""
    data = result.get(_TABLE.get(app)) if app in _TABLE else None
    if not data or not data.get("rows"):
        return None
    return data["columns"], data["rows"]


def series(result, app, columns=None, max_points=None):
    """Chosen columns of a dynamic run, decimated to at most max_points.

    ``columns``: names or glob patterns (``*_speed``); the time column is
    always first.  Keeps the first and last rows.
    """
    t = table(result, app)
    if t is None:
        raise ValueError("a %s job has no time series" % app)
    cols, rows = t
    picked = _pick(cols, columns, _TIME[app])
    idx = [cols.index(c) for c in picked]
    n = len(rows)
    stride = 1
    if max_points and n > max_points:
        stride = math.ceil((n - 1) / (max_points - 1)) if max_points > 1 else n
    keep = list(range(0, n, stride))
    if keep[-1] != n - 1:
        keep.append(n - 1)
    return {"columns": picked, "rows": [[rows[i][j] for j in idx] for i in keep],
            "points": len(keep), "decimated_from": n if stride > 1 else None}


def _pick(cols, patterns, time):
    if not patterns:
        return list(cols)
    out = [time]
    for p in patterns:
        hits = [c for c in cols if c == p or fnmatch.fnmatchcase(c, p)]
        if not hits:
            raise ValueError("no column matches %r; columns are %s"
                             % (p, ", ".join(cols)))
        out += [h for h in hits if h not in out]
    return out


def events_from_xml(workdir):
    """The staged input's dynamic events as ``[{"kind", "start", "end"}]``
    or ``[{"kind", "time"}]``; [] when there are none."""
    out = []
    for xml in sorted(Path(workdir).glob("*.xml")):
        try:
            root = ET.parse(xml).getroot()
        except ET.ParseError:
            continue
        if root.tag != "Configuration":
            continue
        for e in root.findall("Dynamic_simulation/Events/*"):
            g = {c.tag: (c.text or "").strip() for c in e}
            if e.tag == "faultEvent":
                out.append({"kind": "fault", "start": float(g["beginFault"]),
                            "end": float(g["endFault"]),
                            "where": "bus %s" % g.get("faultBranch", "?").split()[0]})
            elif e.tag == "BusFault":
                out.append({"kind": "fault", "start": float(g["begin"]),
                            "end": float(g["end"]), "where": "bus %s" % g.get("bus")})
            elif e.tag == "LineStatus":
                out.append({"kind": "branch %s" % ("in" if g.get("status", "0") == "1" else "out"),
                            "time": float(g["time"]), "where": "line %s" % g.get("line")})
            elif e.tag == "GenStatus":
                out.append({"kind": "generator %s" % ("in" if g.get("status", "0") == "1" else "out"),
                            "time": float(g["time"]),
                            "where": "bus %s id %s" % (g.get("bus"), g.get("id", "1"))})
        break
    return out


# ---- reports -----------------------------------------------------------

def report(result, app, workdir=None):
    """A few lines a reader who did not run the job can act on."""
    fn = {"powerflow": _report_pf, "se": _report_se, "ca": _report_ca,
          "dsf": _report_dyn, "hadrec": _report_dyn, "emt": _report_emt}
    return "\n".join(fn[app](result, workdir))


def _report_pf(r, workdir):
    buses, branches = r.get("buses", []), r.get("branches", [])
    lines = ["Power flow %s in %d iterations; %d buses, %d branches, %d generators."
             % ("converged" if r.get("converged") else "DID NOT CONVERGE",
                r.get("iterations", 0), len(buses), len(branches),
                len(r.get("generators", [])))]
    if not buses:
        return lines
    gen = sum(b["pGen"] for b in buses)
    load = sum(b["pLoad"] for b in buses)
    lines.append("Generation %.1f MW, load %.1f MW, losses %.1f MW."
                 % (gen, load, sum(b.get("pLoss", 0.0) for b in branches)))
    lo, hi = min(buses, key=lambda b: b["voltage"]), max(buses, key=lambda b: b["voltage"])
    lines.append("Voltage from %.4f p.u. at bus %d to %.4f p.u. at bus %d."
                 % (lo["voltage"], lo["busId"], hi["voltage"], hi["busId"]))
    slack = [b for b in buses if b.get("type") == 3]
    if slack:
        lines.append("Slack bus %d supplies %.1f MW, %.1f MVAr."
                     % (slack[0]["busId"], slack[0]["pGen"], slack[0]["qGen"]))
    lines += _violations(buses, branches)
    return lines


def _violations(buses, branches):
    bad = [b for b in buses if not _VMIN <= b["voltage"] <= _VMAX]
    rated = [b for b in branches if b.get("rateA", 0) > 0]
    over = [b for b in rated if b.get("loadingPercent", 0) > _LOAD]
    out = ["Voltage violations (outside %.2f-%.2f p.u.): %s." % (
        _VMIN, _VMAX, ", ".join("bus %d at %.4f" % (b["busId"], b["voltage"])
                                for b in sorted(bad, key=lambda b: b["voltage"])[:10])
        or "none")]
    if rated:
        out.append("Overloads (> %.0f%% of rate A, %d rated branches): %s." % (
            _LOAD, len(rated), ", ".join(
                "%d-%d at %.0f%%" % (b["fromBus"], b["toBus"], b["loadingPercent"])
                for b in sorted(over, key=lambda b: -b["loadingPercent"])[:10])
            or "none"))
    else:
        out.append("No branch ratings in the case, so overloads cannot be judged.")
    return out


def _report_se(r, workdir):
    bad = r.get("bad_measurements", [])
    lines = ["State estimation %s; %d bus states." % (
        "converged" if r.get("converged") else "DID NOT CONVERGE",
        len(r.get("buses", [])))]
    lines.append("Bad measurements flagged: %d%s." % (
        len(bad), ": " + ", ".join(str(m) for m in bad[:10]) if bad else ""))
    return lines


def _report_ca(r, workdir):
    cs = r.get("contingencies", [])
    counts = {}
    for c in cs:
        counts[c["status"]] = counts.get(c["status"], 0) + 1
    lines = ["Contingency analysis: base case %s; %d contingencies, voltage band %.2f-%.2f p.u."
             % ("converged" if r.get("converged") else "DID NOT CONVERGE", len(cs),
                r.get("min_voltage", 0.0), r.get("max_voltage", 0.0))]
    lines.append("By status: %s." % ", ".join("%s %d" % kv for kv in sorted(counts.items())))
    bad = [c["name"] for c in cs if c["status"] != "ok"]
    if bad:
        lines.append("Not ok: %s%s." % (", ".join(bad[:20]),
                                        " ..." if len(bad) > 20 else ""))
    return lines


def _report_emt(r, workdir):
    return ["EMT run finished (exit %s); the EMT application reports status only."
            % r.get("exit_code")]


def dynamic_stats(result, app, after=0.0):
    """Swing statistics of a dynamic run (speed deviations, damping, rotor
    angle spread in degrees), or None without a time series.  ``after``:
    the first disturbance, before which peaks are not counted."""
    t = table(result, app)
    if t is None:
        return None
    cols, rows = t
    time = [row[0] for row in rows]
    out = {"duration": time[-1], "rows": len(rows), "channels": len(cols) - 1,
           "time_step": time[1] - time[0] if len(time) > 1 else 0.0}
    speeds = [c for c in cols if c.endswith("_speed") or c.endswith("_rspeed")]
    angles = [c for c in cols if c.endswith("_PowerAngle") or c.endswith("_rangle")]
    if speeds:
        dev = {c: [abs(row[cols.index(c)] - 1.0) for row in rows] for c in speeds}
        worst = max(speeds, key=lambda c: max(dev[c]))
        out.update(peak_speed_deviation=max(dev[worst]), worst_channel=worst,
                   final_speed_deviation=dev[worst][-1],
                   damping_ratio=_damping(dev[worst], time, after))
    if len(angles) > 1:
        # Each machine's angle relative to its own start: the watch's
        # absolute PowerAngle values carry arbitrary offsets per machine.
        idx = [cols.index(c) for c in angles]
        base = [rows[0][i] for i in idx]
        spread = [_DEG * (max(row[i] - b for i, b in zip(idx, base))
                          - min(row[i] - b for i, b in zip(idx, base)))
                  for row in rows]
        out.update(angle_spread_peak=max(spread), angle_spread_final=spread[-1],
                   loss_of_synchronism=max(spread) > 180.0)
    if "final_speed_deviation" in out:
        out["runaway_speed"] = out["final_speed_deviation"] > 0.02
    return out


def _report_dyn(r, workdir):
    app = "dsf" if "watch" in r else "hadrec"
    events = events_from_xml(workdir) if workdir else []
    first = min([e.get("start", e.get("time", 0.0)) for e in events] or [0.0])
    st = dynamic_stats(r, app, first)
    if st is None:
        return ["Dynamic run finished but recorded no time series (no generator "
                "watch or observations in the case)."]
    lines = ["Dynamic simulation of %.2f s at %.4g s steps, %d rows, %d channels."
             % (st["duration"], st["time_step"], st["rows"], st["channels"])]
    for e in events:
        if "time" in e:
            lines.append("Event: %s at %.3g s (%s)." % (e["kind"], e["time"], e["where"]))
        else:
            lines.append("Event: %s at %s from %.3g to %.3g s (%.0f ms)."
                         % (e["kind"], e["where"], e["start"], e["end"],
                            1000 * (e["end"] - e["start"])))
    verdict = []
    if "peak_speed_deviation" in st:
        lines.append("Peak speed deviation %.3f%% on %s; final %.3f%%."
                     % (100 * st["peak_speed_deviation"], _gen(st["worst_channel"]),
                        100 * st["final_speed_deviation"]))
        zeta = st["damping_ratio"]
        if zeta is not None:
            lines.append("Estimated damping ratio of that swing: %.1f%%." % (100 * zeta))
            verdict.append("poorly damped" if zeta < 0.03 else "damped")
        verdict.append("speed within 1%" if st["peak_speed_deviation"] < 0.01
                       else "large speed excursion")
    if "angle_spread_peak" in st:
        lines.append("Rotor angle spread between generators (each from its own "
                     "start): peak %.1f deg, final %.1f deg." % (
                         st["angle_spread_peak"], st["angle_spread_final"]))
        verdict.append("angles separated beyond 180 deg: loss of synchronism"
                       if st["loss_of_synchronism"] else "angles held together")
    if st.get("runaway_speed"):
        verdict.append("speed still %.1f%% off at the end: machine running away"
                       % (100 * st["final_speed_deviation"]))
    if verdict:
        lines.append("Verdict: %s." % "; ".join(verdict))
    return lines


def _gen(col):
    bus, gid = col.split("_")[:2] if not col.startswith("gen_") else col.split("_")[1:3]
    return "generator %s at bus %s" % (gid, bus)


def _damping(dev, time, after):
    """Damping ratio from successive peaks of |speed - 1| after ``after``,
    or None with fewer than two peaks."""
    peaks = [dev[i] for i in range(1, len(dev) - 1)
             if time[i] > after and dev[i] >= dev[i - 1] and dev[i] > dev[i + 1]
             and dev[i] > 1e-6]
    if len(peaks) < 2 or peaks[-1] <= 0 or peaks[-2] <= 0:
        return None
    delta = math.log(peaks[-2] / peaks[-1])      # log decrement, one half cycle
    if delta <= 0:
        return 0.0
    return delta / math.sqrt(math.pi ** 2 + delta ** 2)


# ---- comparison --------------------------------------------------------

def compare(a, b, app_a, app_b):
    """What changed from job ``a`` to job ``b`` of the same kind."""
    if (app_a in _TABLE) != (app_b in _TABLE):
        raise ValueError("cannot compare a %s job with a %s job" % (app_a, app_b))
    if app_a in _TABLE:
        return _compare_dyn(a, b, app_a, app_b)
    return _compare_static(a, b)


def _compare_static(a, b):
    out = {}
    ba = {x["busId"]: x for x in a.get("buses", [])}
    bb = {x["busId"]: x for x in b.get("buses", [])}
    common = sorted(set(ba) & set(bb))
    if not common:
        raise ValueError("the two jobs share no buses")
    out["buses"] = {"common": len(common), "only_in_a": sorted(set(ba) - set(bb)),
                    "only_in_b": sorted(set(bb) - set(ba))}
    for key, label in (("voltage", "voltage_pu"), ("angle", "angle_deg")):
        diffs = [(bb[i][key] - ba[i][key], i) for i in common if key in ba[i]]
        if diffs:
            d, bus = max(diffs, key=lambda x: abs(x[0]))
            out[label] = {"max_change": d, "at_bus": bus,
                          "rms_change": math.sqrt(sum(x * x for x, _ in diffs) / len(diffs))}
    for key in ("pGen", "pLoad", "qGen", "qLoad"):
        if key in ba[common[0]]:
            out[key + "_total"] = {"a": sum(ba[i][key] for i in common),
                                   "b": sum(bb[i][key] for i in common)}
    ka = {(x["fromBus"], x["toBus"], x["circuitId"].strip()): x for x in a.get("branches", [])}
    kb = {(x["fromBus"], x["toBus"], x["circuitId"].strip()): x for x in b.get("branches", [])}
    out["branches"] = {"only_in_a": ["%d-%d %s" % k for k in sorted(set(ka) - set(kb))],
                       "only_in_b": ["%d-%d %s" % k for k in sorted(set(kb) - set(ka))]}
    flows = sorted(((kb[k]["pFrom"] - ka[k]["pFrom"], k) for k in set(ka) & set(kb)),
                   key=lambda x: -abs(x[0]))[:5]
    out["largest_flow_changes_mw"] = [{"branch": "%d-%d %s" % k, "change": d}
                                      for d, k in flows]
    if "converged" in a or "converged" in b:
        out["converged"] = {"a": a.get("converged"), "b": b.get("converged")}
    return out


def _compare_dyn(a, b, app_a, app_b):
    ta, tb = table(a, app_a), table(b, app_b)
    if ta is None or tb is None:
        raise ValueError("both jobs need a recorded time series")
    (ca, ra), (cb, rb) = ta, tb
    common = [c for c in ca if c in cb and c != ca[0]]
    if not common:
        raise ValueError("the two runs share no channels")
    n = min(len(ra), len(rb))
    if any(abs(ra[i][0] - rb[i][0]) > 1e-9 for i in range(n)):
        raise ValueError("the two runs use different time steps")
    chans = []
    for c in common:
        ia, ib = ca.index(c), cb.index(c)
        diffs = [rb[i][ib] - ra[i][ia] for i in range(n)]
        k = max(range(n), key=lambda i: abs(diffs[i]))
        chans.append({"channel": c, "max_change": diffs[k], "at_time": ra[k][0],
                      "final_change": diffs[n - 1]})
    chans.sort(key=lambda x: -abs(x["max_change"]))
    first = next((ra[i][0] for i in range(n)
                  if any(abs(rb[i][cb.index(c)] - ra[i][ca.index(c)]) > 1e-9
                         for c in common)), None)
    return {"rows_compared": n, "rows": {"a": len(ra), "b": len(rb)},
            "first_difference_at": first, "channels": chans}
