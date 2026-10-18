# -------------------------------------------------------------
# file: python/tests/test_post.py
# -------------------------------------------------------------
# gridpack.mcp.post and gridpack.mcp.plots on hand-made results.
# -------------------------------------------------------------

import math
import re
import sys

import pytest

from gridpack.mcp import plots, post


def _dyn(n=101, step=0.01, fault=(0.1, 0.15), damp=2.0):
    """A watch table: gen 1 swings after the fault and decays; gen 2 is
    flat; angles in radians."""
    cols = ["t", "1_1_speed", "1_1_PowerAngle", "2_1_speed", "2_1_PowerAngle"]
    rows = []
    for i in range(n):
        t = round(i * step, 6)
        x = max(0.0, t - fault[1])
        swing = 0.004 * math.exp(-damp * x) * math.cos(12 * x) if t > fault[1] else 0.0
        rows.append([t, 1.0 + swing, 0.5 + 10 * swing, 1.0, 0.3])
    return {"time_step": step, "final_time": rows[-1][0],
            "watch": {"columns": cols, "rows": rows}}


def _pf(v9=1.0, open_15=False):
    buses = [{"busId": b, "name": "B%d" % b, "type": 3 if b == 1 else 1,
              "voltage": v9 if b == 9 else 1.02, "angle": -float(b),
              "pGen": 300.0 if b == 1 else 0.0, "qGen": 10.0 if b == 1 else 0.0,
              "pLoad": 0.0 if b == 1 else 20.0, "qLoad": 5.0, "shuntMvar": 0.0}
             for b in range(1, 15)]
    branches = [{"fromBus": 1, "toBus": 2, "circuitId": "BL", "pFrom": 150.0,
                 "pLoss": 4.0, "rateA": 100.0, "loadingPercent": 120.0},
                {"fromBus": 1, "toBus": 5, "circuitId": "BL", "pFrom": 75.0,
                 "pLoss": 2.0, "rateA": 0.0, "loadingPercent": 0.0}]
    if open_15:
        branches = branches[:1]
    return {"converged": True, "iterations": 3, "buses": buses,
            "branches": branches, "generators": [{"busId": 1, "genId": "1 "}]}


def test_series_selects_and_decimates():
    r = _dyn()
    s = post.series(r, "dsf", ["*_speed"], max_points=11)
    assert s["columns"] == ["t", "1_1_speed", "2_1_speed"]
    assert (s["points"], s["decimated_from"]) == (11, 101)
    assert [row[0] for row in s["rows"]] == [round(0.1 * k, 6) for k in range(11)]
    full = post.series(r, "dsf")
    assert full["columns"] == r["watch"]["columns"] and full["decimated_from"] is None
    # The last row survives an uneven stride.
    assert post.series(r, "dsf", max_points=7)["rows"][-1][0] == 1.0
    with pytest.raises(ValueError, match="no column matches 'nope'"):
        post.series(r, "dsf", ["nope"])
    with pytest.raises(ValueError, match="has no time series"):
        post.series(_pf(), "powerflow")


def test_report_dynamic_verdicts(tmp_path):
    (tmp_path / "in.xml").write_text(
        "<Configuration><Dynamic_simulation><Events>"
        "<faultEvent><beginFault>0.1</beginFault><endFault>0.15</endFault>"
        "<faultBranch>6 7</faultBranch></faultEvent>"
        "<GenStatus><time>0.5</time><bus>3</bus><id>1</id><status>0</status>"
        "</GenStatus></Events></Dynamic_simulation></Configuration>")
    text = post.report(_dyn(), "dsf", tmp_path)
    assert "1.00 s at 0.01 s steps, 101 rows, 4 channels" in text
    assert "Event: fault at bus 6 from 0.1 to 0.15 s (50 ms)" in text
    assert "Event: generator out at 0.5 s (bus 3 id 1)" in text
    assert re.search(r"Peak speed deviation 0\.3\d\d% on generator 1 at bus 1; final 0\.0\d\d%", text)
    assert "damping ratio" in text and "Verdict: damped; speed within 1%" in text
    assert "angles held together" in text
    # Growing angle separation is called out.
    r = _dyn()
    for i, row in enumerate(r["watch"]["rows"]):
        row[2] = 0.5 + 0.05 * i                       # to 5.5 rad = 315 deg
    assert "loss of synchronism" in post.report(r, "dsf", None)
    assert post.report({"watch": None}, "dsf", None).startswith(
        "Dynamic run finished but recorded no time series")


def test_damping_estimate_matches_the_decay():
    for damp, zeta_true in ((1.5, 0.124), (6.0, 0.447)):
        r = _dyn(n=401, damp=damp)
        dev = [abs(row[1] - 1.0) for row in r["watch"]["rows"]]
        t = [row[0] for row in r["watch"]["rows"]]
        # exp(-d x) cos(12 x): zeta = d / sqrt(d^2 + 12^2)
        assert abs(post._damping(dev, t, 0.15) - zeta_true) < 0.02, damp
    assert post._damping([0.0] * 10, [0.01 * i for i in range(10)], 0.0) is None


def test_report_static_and_others():
    text = post.report(_pf(v9=0.85), "powerflow", None)
    assert "converged in 3 iterations; 14 buses, 2 branches, 1 generators" in text
    assert "Generation 300.0 MW, load 260.0 MW, losses 6.0 MW" in text
    assert "Voltage from 0.8500 p.u. at bus 9" in text
    assert "Slack bus 1 supplies 300.0 MW, 10.0 MVAr" in text
    assert "Voltage violations (outside 0.90-1.10 p.u.): bus 9 at 0.8500" in text
    assert "Overloads (> 100% of rate A, 1 rated branches): 1-2 at 120%" in text
    assert "DID NOT CONVERGE" in post.report(dict(_pf(), converged=False), "powerflow", None)
    ca = {"converged": True, "min_voltage": 0.9, "max_voltage": 1.1,
          "contingencies": [{"name": "L1", "status": "ok"},
                            {"name": "L2", "status": "voltage violation"}]}
    text = post.report(ca, "ca", None)
    assert "2 contingencies" in text and "Not ok: L2" in text
    assert "ok 1, voltage violation 1" in text
    assert "Bad measurements flagged: 2: m1, m2" in post.report(
        {"converged": True, "buses": [1] * 14, "bad_measurements": ["m1", "m2"]},
        "se", None)
    assert "status only" in post.report({"exit_code": 0}, "emt", None)


def test_compare_static_and_dynamic():
    c = post.compare(_pf(), _pf(v9=0.95, open_15=True), "powerflow", "powerflow")
    assert c["voltage_pu"] == {"max_change": pytest.approx(-0.05), "at_bus": 9,
                               "rms_change": pytest.approx(0.05 / math.sqrt(14))}
    assert c["branches"] == {"only_in_a": ["1-5 BL"], "only_in_b": []}
    assert c["pLoad_total"] == {"a": 260.0, "b": 260.0}
    a, b = _dyn(), _dyn(damp=0.5)
    d = post.compare(a, b, "dsf", "dsf")
    assert d["rows_compared"] == 101 and d["first_difference_at"] == 0.16
    assert d["channels"][0]["channel"] in ("1_1_PowerAngle", "1_1_speed")
    assert all(ch["max_change"] == 0 for ch in d["channels"] if ch["channel"].startswith("2_"))
    with pytest.raises(ValueError, match="cannot compare a powerflow job with a dsf"):
        post.compare(_pf(), _dyn(), "powerflow", "dsf")
    with pytest.raises(ValueError, match="different time steps"):
        post.compare(_dyn(), _dyn(step=0.02), "dsf", "dsf")


def test_events_from_xml_ignores_other_xml(tmp_path):
    (tmp_path / "contingencies.xml").write_text("<Contingency_List/>")
    (tmp_path / "in.xml").write_text(
        "<Configuration><Dynamic_simulation><Events><BusFault><begin>1</begin>"
        "<end>1.1</end><bus>7</bus></BusFault><LineStatus><time>2</time>"
        "<line>5 1</line><id>1</id><status>1</status></LineStatus></Events>"
        "</Dynamic_simulation></Configuration>")
    assert post.events_from_xml(tmp_path) == [
        {"kind": "fault", "start": 1.0, "end": 1.1, "where": "bus 7"},
        {"kind": "branch in", "time": 2.0, "where": "line 5 1"}]
    assert post.events_from_xml(tmp_path / "nowhere") == []


def test_plot_without_matplotlib_says_how_to_install(monkeypatch):
    monkeypatch.setitem(sys.modules, "matplotlib", None)
    with pytest.raises(ValueError, match="pip install 'gridpack\\[plot\\]'"):
        plots.render(_dyn(), "dsf")


def test_plot_renders_png():
    pytest.importorskip("matplotlib")
    png = plots.render(_dyn(), "dsf", channels=["*_speed", "1_1_PowerAngle"],
                       events=[{"kind": "fault", "start": 0.1, "end": 0.15},
                               {"kind": "generator out", "time": 0.5}])
    assert png[:8] == b"\x89PNG\r\n\x1a\n" and len(png) > 10000
    assert plots.render(_pf(), "powerflow")[:4] == b"\x89PNG"
    assert plots.render(_pf(), "powerflow", kind="angle")[:4] == b"\x89PNG"
    for kwargs, msg in [({"channels": ["nope"]}, "no channel matches"),
                        ({"kind": "bars"}, 'plots kind "series"')]:
        with pytest.raises(ValueError, match=msg):
            plots.render(_dyn(), "dsf", **kwargs)
    with pytest.raises(ValueError, match="no plot for a ca job"):
        plots.render({}, "ca")
    with pytest.raises(ValueError, match="too many"):
        big = _dyn()
        big["watch"]["columns"] += ["%d_1_speed" % k for k in range(3, 12)]
        for row in big["watch"]["rows"]:
            row += [1.0] * 9
        plots.render(big, "dsf", channels=["*_speed"])
