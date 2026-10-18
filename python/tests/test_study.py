# -------------------------------------------------------------
# file: python/tests/test_study.py
# -------------------------------------------------------------
# gridpack.mcp.study: branch lists and ranking, without GridPACK.
# -------------------------------------------------------------

import pytest

from gridpack.mcp import study
from .conftest import _data_sets


class _Job:
    def __init__(self, state, result=None, error=None):
        self.state, self.result, self.error = state, result, error


class _Runner:
    def __init__(self, jobs):
        self.jobs = jobs

    def get(self, job_id):
        return self.jobs[job_id]


def _watch(spread_deg, speed_dev, n=20):
    import math
    cols = ["t", "1_1_speed", "1_1_angle", "2_1_speed", "2_1_angle"]
    rows = [[0.01 * i, 1.0 + speed_dev * i / n, 0.0, 1.0,
             math.radians(spread_deg) * i / n] for i in range(n + 1)]
    return {"watch": {"columns": cols, "rows": rows}}


def test_n1_branches_lists_and_filters():
    raw = _data_sets() / "raw/kundur-twoarea_v33.raw"
    allb = study.n1_branches(raw)
    assert len(allb) == 12 and (5, 1, "1") in allb       # 8 lines + 4 transformers
    assert study.n1_branches(raw, ["8-7 2", "5-1"]) == [(7, 8, "2"), (5, 1, "1")]
    for only, msg in [(["7-8"], "circuits 1, 2; name one"),
                      (["1-9"], "no branch 1-9"), (["x"], "not \"from-to\"")]:
        with pytest.raises(ValueError, match=msg):
            study.n1_branches(raw, only)


def test_rank_dynamic_n1_most_severe_first():
    st = study.Study("dynamic_n1", "base.xml", {"fault_start": 0.1})
    jobs = {"a": _Job("done", _watch(30, 0.001)),
            "f": _Job("done", _watch(90, 0.001)),      # same spread as c:
            "b": _Job("done", _watch(200, 0.002)),     # speed breaks the tie
            "c": _Job("done", _watch(90, 0.004)),
            "d": _Job("failed", None, "exit 59"),
            "e": _Job("running")}
    for k in jobs:
        st.items.append({"name": "out-" + k, "case": "derived/%s/x.xml" % k,
                         "job_id": k, "branch": k})
    r = study.rank(st, _Runner(jobs))
    assert [row["name"] for row in r["table"]] == ["out-b", "out-c", "out-f",
                                                    "out-a"]
    assert r["table"][0]["unstable"] and not r["table"][1]["unstable"]
    # A machine still 3% fast at the end counts as unstable too.
    jobs["g"] = _Job("done", _watch(10, 0.03))
    st.items.append({"name": "out-g", "case": "c", "job_id": "g", "branch": "g"})
    assert [row["name"] for row in study.rank(st, _Runner(jobs))["table"]][:2] == [
        "out-b", "out-g"]
    del jobs["g"]; st.items.pop()
    assert r["table"][1]["angle_spread_peak_deg"] == pytest.approx(90.0)
    assert (r["finished"], r["pending"]) == (4, 1)
    assert r["failed"] == [{"name": "out-d", "job_id": "d", "error": "exit 59"}]
    assert st.status(_Runner(jobs))["states"] == {"done": 4, "failed": 1,
                                                   "running": 1}


def test_rank_load_sweep_by_factor():
    def pf(factor, v):
        return {"converged": True, "iterations": 4,
                "buses": [{"busId": 1, "type": 3, "voltage": 1.0, "pGen": 100 * factor,
                           "pLoad": 0.0},
                          {"busId": 2, "type": 1, "voltage": v, "pGen": 0.0,
                           "pLoad": 95 * factor}],
                "branches": [{"pLoss": 5 * factor, "rateA": 100.0,
                              "loadingPercent": 90 * factor}]}
    st = study.Study("load_sweep", "pf.xml", {"factors": [1.2, 1.0]})
    jobs = {"x": _Job("done", pf(1.2, 0.88)), "y": _Job("done", pf(1.0, 0.97))}
    st.items = [{"name": "load-1.2", "case": "c1", "job_id": "x", "factor": 1.2},
                {"name": "load-1", "case": "c2", "job_id": "y", "factor": 1.0}]
    t = study.rank(st, _Runner(jobs))["table"]
    assert [r["factor"] for r in t] == [1.0, 1.2]
    assert t[1] == {"factor": 1.2, "converged": True, "iterations": 4,
                    "load_mw": 114.0, "generation_mw": 120.0, "losses_mw": 6.0,
                    "min_voltage": 0.88, "min_voltage_bus": 2, "slack_mw": 120.0,
                    "voltage_violations": 1, "max_loading_percent": 108.0,
                    "name": "load-1.2", "case": "c1", "job_id": "x"}
