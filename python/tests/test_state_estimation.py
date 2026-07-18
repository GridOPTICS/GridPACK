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
            se = StateEstimation(s, "input.xml", suppress_output=True)
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
            se = StateEstimation(s, "input.xml", suppress_output=True)
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
