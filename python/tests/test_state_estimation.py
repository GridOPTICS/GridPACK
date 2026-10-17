# -------------------------------------------------------------
# file: python/tests/test_state_estimation.py
# -------------------------------------------------------------
# State estimation integration tests: IEEE 14 bus + measurement file.
# -------------------------------------------------------------

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
