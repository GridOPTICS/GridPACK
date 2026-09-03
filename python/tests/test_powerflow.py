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


# -------------------------------------------------------------
# Convergence reporting
# -------------------------------------------------------------
# Divergence is manufactured by capping maxIteration rather than shipping
# an ill-conditioned case; what is under test is the reporting path.


def _capped_case(tests_data_dir, tmp_path, max_iteration=1):
    """Copy input_14.xml into tmp_path with maxIteration capped."""
    import shutil

    shutil.copy(tests_data_dir / "IEEE14.raw", tmp_path / "IEEE14.raw")
    xml = (tests_data_dir / "input_14.xml").read_text()
    assert "<maxIteration>50</maxIteration>" in xml, (
        "input_14.xml no longer has maxIteration=50; update this helper"
    )
    xml = xml.replace(
        "<maxIteration>50</maxIteration>",
        f"<maxIteration>{max_iteration}</maxIteration>",
    )
    (tmp_path / "capped_14.xml").write_text(xml)
    return tmp_path


@pytest.mark.integration
def test_converged_result_reports_convergence(tests_data_dir):
    """A converged solve exposes summary, iterations and mismatch."""
    r = run_inline(
        """
        from gridpack import Session, PowerFlow
        with Session() as s:
            pf = PowerFlow(s, "input_14.xml", suppress_output=True)
            result = pf.solve()
            print(f"CONVERGED={result.converged}")
            print(f"ITERATIONS={result.iterations}")
            print(f"HAS_SUMMARY={result.convergence is not None}")
            print(f"PER_ITER={len(result.convergence.perIteration)}")
            print(f"WORST_Q_BUS={result.mismatch.maxQBus}")
        """,
        cwd=tests_data_dir,
    )
    assert r.returncode == 0, f"stdout:\n{r.stdout}\nstderr:\n{r.stderr}"
    assert "CONVERGED=True" in r.stdout
    assert "ITERATIONS=2" in r.stdout
    assert "HAS_SUMMARY=True" in r.stdout
    assert "PER_ITER=2" in r.stdout


@pytest.mark.integration
def test_diverged_raises_by_default(tests_data_dir, tmp_path):
    """strict defaults to True, so a capped solve raises."""
    d = _capped_case(tests_data_dir, tmp_path)
    r = run_inline(
        """
        from gridpack import Session, PowerFlow, PowerFlowDiverged, GridPACKError
        with Session() as s:
            pf = PowerFlow(s, "capped_14.xml", suppress_output=True)
            try:
                pf.solve()
                print("RAISED=no")
            except PowerFlowDiverged as e:
                print(f"RAISED=yes")
                print(f"IS_GRIDPACK_ERROR={isinstance(e, GridPACKError)}")
                print(f"ITERATIONS={e.iterations}")
                print(f"TOLERANCE={e.tolerance}")
                print(f"MAX_ITERATION={e.max_iteration}")
                print(f"HAS_MISMATCH={e.mismatch is not None}")
                # message must carry numbers, not just "diverged"
                m = str(e)
                print(f"MSG_HAS_TOL={'tolerance' in m}")
                print(f"MSG_HAS_BUS={'bus' in m}")
        """,
        cwd=d,
    )
    assert r.returncode == 0, f"stdout:\n{r.stdout}\nstderr:\n{r.stderr}"
    assert "RAISED=yes" in r.stdout
    assert "IS_GRIDPACK_ERROR=True" in r.stdout
    assert "ITERATIONS=1" in r.stdout
    assert "MAX_ITERATION=1" in r.stdout
    assert "HAS_MISMATCH=True" in r.stdout
    assert "MSG_HAS_TOL=True" in r.stdout
    assert "MSG_HAS_BUS=True" in r.stdout


@pytest.mark.integration
def test_diverged_strict_false_returns_result(tests_data_dir, tmp_path):
    """strict=False hands back the unconverged result for inspection."""
    d = _capped_case(tests_data_dir, tmp_path)
    r = run_inline(
        """
        from gridpack import Session, PowerFlow
        with Session() as s:
            pf = PowerFlow(s, "capped_14.xml", suppress_output=True)
            result = pf.solve(strict=False)
            print(f"CONVERGED={result.converged}")
            print(f"ITERATIONS={result.iterations}")
            print(f"WORST_Q_BUS={result.mismatch.maxQBus}")
            print(f"REPR_SAYS_DIVERGED={'DIVERGED' in repr(result)}")
        """,
        cwd=d,
    )
    assert r.returncode == 0, f"stdout:\n{r.stdout}\nstderr:\n{r.stderr}"
    assert "CONVERGED=False" in r.stdout
    assert "ITERATIONS=1" in r.stdout
    assert "REPR_SAYS_DIVERGED=True" in r.stdout


@pytest.mark.integration
def test_nonlinear_path_reports_via_solver_return(tests_data_dir):
    """nl_solve populates no summary, so the verdict is its return value.

    Reading the summary here would report converged=False for a run that
    converged.
    """
    r = run_inline(
        """
        from gridpack import Session, PowerFlow
        with Session() as s:
            pf = PowerFlow(s, "input_14.xml", suppress_output=True)
            result = pf.solve(nonlinear=True)      # strict=True must not fire
            print(f"CONVERGED={result.converged}")
            print(f"SUMMARY_IS_NONE={result.convergence is None}")
            print(f"ITERATIONS_IS_NONE={result.iterations is None}")
            print(f"MISMATCH_IS_NONE={result.mismatch is None}")
        """,
        cwd=tests_data_dir,
    )
    assert r.returncode == 0, f"stdout:\n{r.stdout}\nstderr:\n{r.stderr}"
    assert "CONVERGED=True" in r.stdout
    assert "SUMMARY_IS_NONE=True" in r.stdout
    assert "ITERATIONS_IS_NONE=True" in r.stdout
    assert "MISMATCH_IS_NONE=True" in r.stdout


@pytest.mark.integration
@pytest.mark.mpi
def test_diverged_raises_on_every_rank(tests_data_dir, tmp_path, require_mpiexec):
    """Every rank must raise, or MPI code deadlocks in the handler.

    The verdict is globally reduced, so diagnostics must match per rank.
    """
    d = _capped_case(tests_data_dir, tmp_path)
    r = run_inline(
        """
        from gridpack import Session, PowerFlow, PowerFlowDiverged
        with Session() as s:
            pf = PowerFlow(s, "capped_14.xml", suppress_output=True)
            try:
                pf.solve()
                print(f"RANK{s.rank}=no-raise")
            except PowerFlowDiverged as e:
                print(f"RANK{s.rank}=raised tol={e.final_tolerance:.12e} "
                      f"iters={e.iterations} qbus={e.mismatch.maxQBus}")
        """,
        cwd=d,
        mpi_np=2,
    )
    assert r.returncode == 0, f"stdout:\n{r.stdout}\nstderr:\n{r.stderr}"
    lines = {}
    for line in r.stdout.splitlines():
        line = line.strip()
        for rank in (0, 1):
            if line.startswith(f"RANK{rank}="):
                lines[rank] = line[len(f"RANK{rank}="):]
    assert set(lines) == {0, 1}, f"missing a rank in:\n{r.stdout}"
    assert "no-raise" not in lines[0] and "no-raise" not in lines[1], (
        f"a rank did not raise -- collective hazard:\n{r.stdout}"
    )
    assert lines[0] == lines[1], (
        f"ranks disagree on the diagnostic:\n  0: {lines[0]}\n  1: {lines[1]}"
    )


@pytest.mark.parametrize(
    "raw,expected",
    [
        ("true", True), ("True", True), ("TRUE", True), ("yes", True),
        ("on", True), ("1", True),
        ("false", False), ("False", False), ("no", False), ("off", False),
        ("0", False), ("", False), (None, False),
        ("  true  ", True), ("  false ", False),
        ("nonsense", False),          # unrecognized -> default
        (True, True), (False, False),
    ],
)
def test_xml_bool_coercion(raw, expected):
    """XML scalars are strings, so bool("false") would be True."""
    from gridpack.powerflow import _xml_bool
    assert _xml_bool(raw) is expected


def test_xml_number_coercion():
    from gridpack.powerflow import _xml_number
    assert _xml_number("1.0e-6", float) == 1.0e-6
    assert _xml_number(" 50 ", int) == 50
    assert _xml_number(None, float) is None
    assert _xml_number("abc", float) is None
    assert _xml_number("abc", float, default=2.0) == 2.0
