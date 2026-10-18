# -------------------------------------------------------------
# file: python/tests/test_measurements.py
# -------------------------------------------------------------
# gridpack.mcp.measurements on a hand-made power flow solution.
# -------------------------------------------------------------

import xml.etree.ElementTree as ET

import pytest

from gridpack.mcp import measurements as M
from .conftest import _data_sets


def _pf():
    buses = [{"busId": b, "voltage": 1.0 + 0.01 * b, "angle": -2.0 * b,
              "pInjection": 10.0 * b, "qInjection": -5.0 * b} for b in (1, 2, 3)]
    branches = [{"fromBus": 1, "toBus": 2, "circuitId": "BL", "pFrom": 50.0,
                 "qFrom": 5.0, "pTo": -49.0, "qTo": -4.0},
                {"fromBus": 2, "toBus": 3, "circuitId": "1 ", "pFrom": 20.0,
                 "qFrom": 2.0, "pTo": -19.5, "qTo": -1.5}]
    return {"converged": True, "buses": buses, "branches": branches}


def test_exact_meters_without_noise():
    meters, summary = M.synthesize(_pf(), 100.0, noise_scale=0.0)
    assert summary["by_type"] == {"VM": 3, "PI": 3, "QI": 3, "PIJ": 2, "QIJ": 2}
    by = {M._key(m): m for m in meters}
    assert by[("VM", 2)]["value"] == 1.02 and by[("PI", 2)]["value"] == 0.2
    assert by[("PIJ", 1, 2, "BL")]["value"] == 0.5 and by[("QIJ", 2, 3, "1")]["value"] == 0.02
    # Exact values still carry the nominal deviation, the estimator's weight.
    assert all(m["deviation"] == M.DEFAULT_STD[m["type"]] for m in meters)
    va, _ = M.synthesize(_pf(), 100.0, types=["VA", "PJI"], noise_scale=0.0)
    assert [m["value"] for m in va] == [-2.0, -4.0, -6.0, -0.49, -0.195]


def test_noise_is_seeded_and_scaled():
    a, _ = M.synthesize(_pf(), 100.0, seed=7)
    b, _ = M.synthesize(_pf(), 100.0, seed=7)
    c, _ = M.synthesize(_pf(), 100.0, seed=8)
    assert [m["value"] for m in a] == [m["value"] for m in b]
    assert [m["value"] for m in a] != [m["value"] for m in c]
    assert all(m["deviation"] == M.DEFAULT_STD[m["type"]] for m in a)
    big, _ = M.synthesize(_pf(), 100.0, seed=7, std={"VM": 0.1}, noise_scale=2.0)
    vm = [m for m in big if m["type"] == "VM"]
    assert all(m["deviation"] == 0.2 for m in vm)
    assert max(abs(m["value"] - m["truth"]) for m in vm) > 0.05


def test_coverage_and_planted_bad_data():
    bad = [{"type": "VM", "bus": 3, "error": 0.2},
           {"type": "PIJ", "from": 2, "to": 3, "ckt": "1 ", "error": -0.5}]
    meters, summary = M.synthesize(_pf(), 100.0, coverage=0.5, seed=1,
                                   noise_scale=0.0, bad_data=bad)
    assert summary["meters"] == round(0.5 * 13) == 6     # planted ones included
    by = {M._key(m): m for m in meters}
    # Planted meters survive any coverage, whatever the seed draws.
    for seed in range(10):
        few, _ = M.synthesize(_pf(), 100.0, coverage=0.25, seed=seed, bad_data=bad)
        assert len(few) == 3 and {("VM", 3), ("PIJ", 2, 3, "1")} <= {M._key(m) for m in few}
    assert by[("VM", 3)]["value"] == pytest.approx(1.23)         # kept and shifted
    assert by[("PIJ", 2, 3, "1")]["value"] == pytest.approx(-0.3)
    assert [p["type"] for p in summary["planted"]] == ["VM", "PIJ"]
    assert summary["planted"][1]["planted_error"] == -0.5
    rnd, s2 = M.synthesize(_pf(), 100.0, random_bad=2, bad_sigma=10, seed=3)
    assert len(s2["planted"]) == 2
    assert all(abs(p["planted_error"]) == 10 * M.DEFAULT_STD[p["type"]]
               for p in s2["planted"])


@pytest.mark.parametrize("kwargs, msg", [
    ({"types": ["VM", "ZZ"]}, "unknown measurement type ZZ"),
    ({"coverage": 0}, "coverage must be"),
    ({"coverage": 1.5}, "coverage must be"),
    ({"bad_data": [{"type": "VM", "bus": 9, "error": 1}]}, "no VM meter"),
    ({"bad_data": [{"type": "VA", "bus": 1, "error": 1}]}, "no VA meter"),
    ({"bad_data": [{"type": "PIJ", "from": 1, "error": 1}]}, "needs from, to and ckt"),
    ({"bad_data": [{"type": "VM", "error": 1}]}, "needs bus"),
    ({"bad_data": [{"type": "IIJ", "from": 1, "to": 2, "ckt": "BL", "error": 1}]},
     "unknown type"),
    ({"bad_data": [{"bus": 1}]}, "must have type and error"),
    ({"random_bad": 99}, "exceeds the 13 meters"),
])
def test_rejections(kwargs, msg):
    with pytest.raises(ValueError, match=msg):
        M.synthesize(_pf(), 100.0, **kwargs)


def test_measurement_file_and_se_xml(tmp_path):
    meters, _ = M.synthesize(_pf(), 100.0, types=["VM", "PIJ"], noise_scale=0.0)
    M.write_measurements(tmp_path / "m.xml", meters)
    root = ET.parse(tmp_path / "m.xml").getroot()
    recs = [{c.tag: c.text for c in m} for m in root.findall("Measurement")]
    assert recs[0] == {"Type": "VM", "Bus": "1", "Value": "1.010000",
                       "Deviation": "0.004"}
    assert recs[3] == {"Type": "PIJ", "FromBus": "1", "ToBus": "2", "CKT": "BL",
                       "Value": "0.500000", "Deviation": "0.008"}
    xml = M.se_xml("/x/net.raw", 33, "m.xml")
    se = ET.fromstring(xml).find("State_estimation")
    assert se.find("networkConfiguration_v33").text == "/x/net.raw"
    assert se.find("measurementList").text == "m.xml"
    assert "<networkConfiguration>" in M.se_xml("n.raw", 23, "m.xml")


def test_raw_sbase_reads_every_header_style():
    raw = _data_sets() / "raw"
    assert M.raw_sbase(raw / "IEEE14.raw") == 100.0                 # v23, spaces
    assert M.raw_sbase(raw / "kundur-twoarea_v33.raw") == 100.0     # commas
    assert M.raw_sbase(raw / "240busWECC_2018_PSS_PQLoad_fixedshunt_noremotebus_yuan_v36.raw") == 100.0


def test_accuracy_and_report_lines():
    truth = _pf()["buses"]
    se = {"buses": [{"busId": 1, "voltage": 1.011, "angle": -2.0},
                    {"busId": 2, "voltage": 1.02, "angle": -3.9},
                    {"busId": 3, "voltage": 1.03, "angle": -6.0}]}
    planted = [{"type": "VM", "bus": 3}, {"type": "PIJ", "from": 1, "to": 2, "ckt": "BL"}]
    flagged = [{"type": "VM", "busId": 3, "value": 1.2, "deviation": 0.004},
               {"type": "PI", "busId": 1, "value": 0.1, "deviation": 0.01}]
    acc = M.accuracy(se, truth, planted, flagged)
    assert acc["voltage_pu"]["max_error"] == pytest.approx(0.001)
    assert acc["voltage_pu"]["at_bus"] == 1
    assert acc["angle_deg"]["max_error"] == pytest.approx(0.1) and acc["angle_deg"]["at_bus"] == 2
    assert acc["bad_data"] == {"planted": 2, "flagged": 2, "caught": ["VM 3"],
                               "missed": ["PIJ 1 2 BL"], "false_alarms": ["PI 1"]}
    lines = M.estimate_lines(acc)
    assert lines[0].startswith("Against the power flow truth: RMS voltage error 0.00058")
    assert lines[1] == ("Bad data: 2 planted, 2 flagged; caught VM 3; missed "
                        "PIJ 1 2 BL; false alarms PI 1.")
