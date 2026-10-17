# -------------------------------------------------------------
# file: python/tests/test_state_estimation.py
# -------------------------------------------------------------
# State estimation integration tests: IEEE 14 bus + measurement file.
# -------------------------------------------------------------

import json

import pytest

from .conftest import run_inline


@pytest.mark.integration
def test_state_estimation_ieee14_serial(se_data_dir):
    r = run_inline(
        """
        from gridpack import Session, StateEstimation
        with Session() as s:
            se = StateEstimation(s, "input_14.xml", suppress_output=True)
            print(f"N_MEAS={len(se.measurements)}")
            result = se.solve()
            if s.rank == 0:
                print(f"CONVERGED={int(result.has_converged())}")
                print("OK")
        """,
        cwd=se_data_dir,
        timeout=180,
    )
    assert r.returncode == 0, (
        f"stdout:\n{r.stdout}\nstderr:\n{r.stderr}"
    )
    assert "OK" in r.stdout
    assert "CONVERGED=1" in r.stdout


@pytest.mark.integration
@pytest.mark.mpi
def test_state_estimation_ieee14_mpi2(se_data_dir, require_mpiexec):
    r = run_inline(
        """
        from gridpack import Session, StateEstimation
        with Session() as s:
            se = StateEstimation(s, "input_14.xml", suppress_output=True)
            result = se.solve()
            if s.rank == 0:
                print(f"CONVERGED={int(result.has_converged())}")
                print("OK")
        """,
        cwd=se_data_dir,
        mpi_np=2,
        timeout=180,
    )
    assert r.returncode == 0, (
        f"stdout:\n{r.stdout}\nstderr:\n{r.stderr}"
    )
    assert "OK" in r.stdout
    assert "CONVERGED=1" in r.stdout


# -------------------------------------------------------------
# Measurement files
# -------------------------------------------------------------
# Measurements used to load through Configuration, which merges every
# file into one process-wide tree and searches only the first: a second
# StateEstimation in the same process silently reused the first file.

def _meas_file(path, vm4):
    path.write_text(
        "<Measurements>"
        "<Measurement><Type>VM</Type><Bus>4</Bus>"
        "<Value>%s</Value><Deviation>0.005</Deviation></Measurement>"
        "<Measurement><Type>PIJ</Type><FromBus>2</FromBus><ToBus>3</ToBus>"
        "<CKT>1</CKT><Value>0.73</Value><Deviation>0.01</Deviation>"
        "</Measurement></Measurements>" % vm4)


def test_read_measurements_matches_the_cpp_fields(tmp_path):
    from gridpack.state_estimation import read_measurements
    _meas_file(tmp_path / "m.xml", "1.0177")
    vm, pij = read_measurements(str(tmp_path / "m.xml"))
    assert (vm.type(), vm.busid, vm.value, vm.deviation) == ("VM", 4, 1.0177, 0.005)
    # CKT is padded to the 2-character PSS/E width, as getMeasurements does.
    assert (pij.type(), pij.fbusid, pij.tbusid, pij.ckt()) == ("PIJ", 2, 3, "1 ")


def test_read_measurements_rejects_a_mismatched_tag(tmp_path):
    """Boost's reader accepts <Values>...</Value> and drops the value."""
    from gridpack.state_estimation import read_measurements
    p = tmp_path / "bad.xml"
    p.write_text("<Measurements><Measurement><Type>VM</Type><Bus>4</Bus>"
                 "<Values>1.0</Value><Deviation>0.005</Deviation>"
                 "</Measurement></Measurements>")
    with pytest.raises(ValueError, match="malformed measurement file"):
        read_measurements(str(p))


@pytest.mark.integration
def test_second_measurement_file_is_not_ignored(se_data_dir, tmp_path):
    import shutil
    for f in se_data_dir.iterdir():
        shutil.copy(f, tmp_path / f.name)
    _meas_file(tmp_path / "a.xml", "1.0177")
    _meas_file(tmp_path / "b.xml", "1.2000")
    r = run_inline("""
        from gridpack import Session, StateEstimation
        with Session() as s:
            for f in ("a.xml", "b.xml"):
                with StateEstimation(s, "input_14.xml", measurement_file=f,
                                     suppress_output=True) as se:
                    print("VM4=%s" % [m.value for m in se.measurements
                                      if m.type() == "VM"])
    """, cwd=tmp_path, timeout=180)
    assert r.returncode == 0, r.stderr[-2000:]
    assert "VM4=[1.0177]" in r.stdout and "VM4=[1.2]" in r.stdout, r.stdout


@pytest.mark.integration
def test_measurements_argument_skips_the_file(se_data_dir):
    r = run_inline("""
        from gridpack import Session, StateEstimation
        from gridpack.state_estimation import read_measurements
        with Session() as s:
            meters = read_measurements("IEEE14_meas.xml")
            with StateEstimation(s, "input_14.xml", measurements=meters,
                                 suppress_output=True) as se:
                n = len(se.measurements)
                ok = se.solve().converged
            print("N=%d CONVERGED=%d" % (n, ok))
    """, cwd=se_data_dir, timeout=180)
    assert r.returncode == 0, r.stderr[-2000:]
    assert "N=66 CONVERGED=1" in r.stdout, r.stdout


# PF truth -> exact meters -> SE; BAD's VM meter reads 12 sigma high.
_PF_TO_SE = """
import json
from gridpack import Session, PowerFlow, StateEstimation
from gridpack.state_estimation import Measurement

def meter(kind, value, **ids):
    m = Measurement(); m.type(kind)
    for k, v in ids.items():
        m.ckt(v) if k == "ckt" else setattr(m, k, v)
    m.value, m.deviation = value, 0.004 if kind == "VM" else 0.01
    return m

with Session() as s:
    with PowerFlow(s, "input_14.xml", suppress_output=True) as pf:
        r = pf.solve(); buses, branches = r.buses(), r.branches()
    meters = []
    for b in buses:
        n = b["busId"]
        meters += [meter("VM", b["voltage"] + (0.05 if n == BAD else 0), busid=n),
                   meter("PI", b["pInjection"] / 100, busid=n),
                   meter("QI", b["qInjection"] / 100, busid=n)]
    for br in branches:
        ids = dict(fbusid=br["fromBus"], tbusid=br["toBus"],
                   ckt=br["circuitId"].strip().ljust(2)[:2])
        meters += [meter("PIJ", br["pFrom"] / 100, **ids),
                   meter("QIJ", br["qFrom"] / 100, **ids)]
    with StateEstimation(s, "input_14.xml", measurements=meters,
                         suppress_output=True) as se:
        res = se.solve()
        est, bad = res.buses(), res.bad_measurements()
    if s.rank == 0:
        print("JSON=" + json.dumps({"truth": buses, "est": est, "bad": bad}))
"""


def _pf_to_se(tmp_path, bad_bus=0, mpi_np=0):
    from gridpack import Case
    case = Case.from_data_sets("input/powerflow/input_14.xml")
    case.include("input/se/input_14.xml", "State_estimation")
    case.override("se", tolerance=1e-8)
    case.set("State_estimation.outputFile", "se_report.txt")
    case.write(tmp_path)
    r = run_inline("BAD = %d\n" % bad_bus + _PF_TO_SE, cwd=tmp_path,
                   mpi_np=mpi_np)
    assert r.returncode == 0, f"stdout:\n{r.stdout}\nstderr:\n{r.stderr}"
    line = next(l for l in r.stdout.splitlines() if l.startswith("JSON="))
    return json.loads(line[5:])


@pytest.mark.integration
def test_estimate_recovers_the_power_flow_state(tmp_path):
    out = _pf_to_se(tmp_path)
    truth = {b["busId"]: b for b in out["truth"]}
    assert [b["busId"] for b in out["est"]] == sorted(truth)
    for b in out["est"]:
        assert abs(b["voltage"] - truth[b["busId"]]["voltage"]) < 1e-5
        assert abs(b["angle"] - truth[b["busId"]]["angle"]) < 1e-4
    assert out["bad"] == []


@pytest.mark.integration
def test_bad_measurements_names_the_corrupted_meter(tmp_path):
    bad = _pf_to_se(tmp_path, bad_bus=12)["bad"]
    assert [(m["type"], m.get("busId")) for m in bad] == [("VM", 12)]


@pytest.mark.integration
@pytest.mark.mpi
def test_estimate_is_the_same_on_two_ranks(tmp_path, require_mpiexec):
    one = _pf_to_se(tmp_path / "np1", bad_bus=12)
    two = _pf_to_se(tmp_path / "np2", bad_bus=12, mpi_np=2)
    assert two["est"] == one["est"] and two["bad"] == one["bad"]


def test_report_parser_reads_bus_and_branch_rows(tmp_path):
    from gridpack.state_estimation import _parse_report
    # Row formats from SEAppModule::write.
    report = tmp_path / "r.txt"
    report.write_text(
        "BAD DATA DETECTION RESULTS\n"
        "Index  Type  Bus/From  To      Value        Deviation    NormResidual\n"
        "  %3d   %-4s  %-8d  -       %-12.6f  %-12.6f  N/A\n"
        % (9, "VM", 12, 1.105, 0.004)
        + "  %3d   %-4s  %-8d  %-7d  %-12.6f  %-12.6f  N/A\n"
        % (40, "PIJ", 2, 3, 1.23, 0.008)
        + "STATE ESTIMATION OUTPUTS - COMPLETE SYSTEM STATE\n"
        "%-10d      %-20.6f      %-20.6f\n" % (2, -4.5, 1.04)
        + "%-10d      %-20.6f      %-20.6f\n" % (1, 0.0, 1.06)
        + "BRANCH POWER FLOW RESULTS - COMPLETE (P.U.)\n"
        "1         2        BL               -4.072177          -17.367337\n")
    assert _parse_report(report) == {
        "buses": [{"busId": 1, "voltage": 1.06, "angle": 0.0},
                  {"busId": 2, "voltage": 1.04, "angle": -4.5}],
        "bad": [{"type": "VM", "busId": 12, "value": 1.105, "deviation": 0.004},
                {"type": "PIJ", "fromBus": 2, "toBus": 3, "value": 1.23,
                 "deviation": 0.008}],
    }
