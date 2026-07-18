# -------------------------------------------------------------
# file: python/tests/test_powerflow.py
# -------------------------------------------------------------
# PowerFlow integration tests: IEEE 14 bus + reference solution.
# Each test drives a subprocess so we get a fresh gridpack.Session
# (Session is a process singleton).
# -------------------------------------------------------------

import pytest

from .conftest import run_inline


BUS5_VMAG_REF = 1.01951386
BUS5_ANGLE_REF = -8.773854
TOL = 1e-5


@pytest.mark.integration
def test_powerflow_ieee14_serial(tests_data_dir):
    r = run_inline(
        """
        from gridpack import Session, PowerFlow
        with Session() as s:
            pf = PowerFlow(s, "input_14.xml", suppress_output=True)
            result = pf.solve()
            vmag, angle = result.get_bus_solution(5)
            print(f"BUS5_VMAG={vmag}")
            print(f"BUS5_ANGLE={angle}")
        """,
        cwd=tests_data_dir,
    )
    assert r.returncode == 0, (
        f"stdout:\n{r.stdout}\nstderr:\n{r.stderr}"
    )
    vmag = _parse_float(r.stdout, "BUS5_VMAG=")
    angle = _parse_float(r.stdout, "BUS5_ANGLE=")
    assert abs(vmag - BUS5_VMAG_REF) < TOL
    assert abs(angle - BUS5_ANGLE_REF) < TOL


@pytest.mark.integration
@pytest.mark.mpi
def test_powerflow_ieee14_mpi2(tests_data_dir, require_mpiexec):
    r = run_inline(
        """
        from gridpack import Session, PowerFlow
        with Session() as s:
            pf = PowerFlow(s, "input_14.xml", suppress_output=True)
            result = pf.solve()
            if s.rank == 0:
                vmag, angle = result.get_bus_solution(5)
                print(f"BUS5_VMAG={vmag}")
                print(f"BUS5_ANGLE={angle}")
        """,
        cwd=tests_data_dir,
        mpi_np=2,
    )
    assert r.returncode == 0, (
        f"stdout:\n{r.stdout}\nstderr:\n{r.stderr}"
    )
    vmag = _parse_float(r.stdout, "BUS5_VMAG=")
    angle = _parse_float(r.stdout, "BUS5_ANGLE=")
    assert abs(vmag - BUS5_VMAG_REF) < TOL
    assert abs(angle - BUS5_ANGLE_REF) < TOL


def _parse_float(text: str, prefix: str) -> float:
    for line in text.splitlines():
        line = line.strip()
        if line.startswith(prefix):
            return float(line[len(prefix):])
    raise AssertionError(f"'{prefix}' not found in output:\n{text}")
