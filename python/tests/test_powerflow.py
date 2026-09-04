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


# -------------------------------------------------------------
# Gathered network tables
# -------------------------------------------------------------
# buses()/branches()/generators() allgather and sort, so their output must
# not depend on the rank count.  The tests below pin that down, and check
# the violation criteria against the C++ checks as an oracle.


def _rated_case(tests_data_dir, tmp_path, rated_raw, scale=None):
    """Build a case in tmp_path from the rated IEEE-14.

    ``scale`` multiplies branch RATEA/B/C; the base case peaks at 84%,
    so scaling down is how an overload is manufactured.
    """
    import re

    raw = rated_raw.read_text().splitlines(True)
    if scale is not None:
        out, in_branch = [], False
        for ln in raw:
            if "BEGIN BRANCH DATA" in ln:
                in_branch = True
            elif "END OF BRANCH DATA" in ln:
                in_branch = False
            elif in_branch and "," in ln:
                f = ln.split(",")
                if len(f) > 8:
                    for i in (6, 7, 8):
                        f[i] = "%8.2f" % (float(f[i]) * scale)
                    ln = ",".join(f)
            out.append(ln)
        raw = out
    (tmp_path / "case.raw").write_text("".join(raw))

    xml = (tests_data_dir / "input_14.xml").read_text()
    xml, n = re.subn(
        r"<networkConfiguration>[^<]*</networkConfiguration>",
        "<networkConfiguration> case.raw </networkConfiguration>",
        xml,
    )
    assert n == 1, "input_14.xml networkConfiguration no longer matches"
    (tmp_path / "case.xml").write_text(xml)
    return tmp_path


_TABLE_DUMP = """
    import json
    from gridpack import Session, PowerFlow
    with Session() as s:
        pf = PowerFlow(s, "{xml}", suppress_output=True)
        r = pf.solve()
        buses, branches, gens = r.buses(), r.branches(), r.generators()
        if s.rank == 0:
            print("COUNTS=%d,%d,%d" % (len(buses), len(branches), len(gens)))
            print("BUSES=" + json.dumps(buses, sort_keys=True))
            print("BRANCHES=" + json.dumps(branches, sort_keys=True))
            print("GENS=" + json.dumps(gens, sort_keys=True))
"""


@pytest.mark.integration
def test_gathered_tables_serial(tests_data_dir):
    """Every bus, branch and generator comes back with all its columns."""
    r = run_inline(_TABLE_DUMP.format(xml="input_14.xml"), cwd=tests_data_dir)
    assert r.returncode == 0, f"stdout:\n{r.stdout}\nstderr:\n{r.stderr}"
    assert "COUNTS=14,20,5" in r.stdout

    import json
    buses = json.loads(_extract(r.stdout, "BUSES="))
    assert {b["busId"] for b in buses} == set(range(1, 15))
    assert [b["busId"] for b in buses] == sorted(b["busId"] for b in buses)
    # Full column set, not just the legacy vmag/vangle view.
    assert set(buses[0]) == {
        "busId", "type", "area", "zone", "baseKV", "voltage", "angle",
        "pInjection", "qInjection", "pLoad", "qLoad", "pGen", "qGen",
        "shuntMvar",
    }
    branches = json.loads(_extract(r.stdout, "BRANCHES="))
    assert set(branches[0]) == {
        "fromBus", "toBus", "circuitId", "pFrom", "qFrom", "pTo", "qTo",
        "pLoss", "qLoss", "mvaFrom", "mvaTo", "rateA", "loadingPercent",
    }


# Parallel linear algebra reassociates, so values move ~1e-14; only row
# identity and ordering are exact.
_RANK_INVARIANT_TOL = 1e-9


@pytest.mark.integration
@pytest.mark.mpi
def test_gathered_tables_rank_invariant(tests_data_dir, require_mpiexec):
    """A correct gather loses and duplicates nothing across rank counts."""
    import json

    body = _TABLE_DUMP.format(xml="input_14.xml")
    serial = run_inline(body, cwd=tests_data_dir)
    par = run_inline(body, cwd=tests_data_dir, mpi_np=2)
    assert serial.returncode == 0, serial.stderr
    assert par.returncode == 0, par.stderr

    assert _extract(serial.stdout, "COUNTS=") == _extract(par.stdout, "COUNTS=")
    for key in ("BUSES=", "BRANCHES=", "GENS="):
        a = json.loads(_extract(serial.stdout, key))
        b = json.loads(_extract(par.stdout, key))
        assert len(a) == len(b), f"{key} row count differs"
        for i, (x, y) in enumerate(zip(a, b)):
            assert set(x) == set(y), f"{key} row {i} column set differs"
            for k in x:
                if isinstance(x[k], float):
                    assert abs(x[k] - y[k]) < _RANK_INVARIANT_TOL, (
                        f"{key} row {i} field {k}: {x[k]} vs {y[k]}"
                    )
                else:
                    # IDs, ordering and circuit tags must match exactly.
                    assert x[k] == y[k], (
                        f"{key} row {i} field {k}: {x[k]!r} vs {y[k]!r}"
                    )


def _extract(text: str, prefix: str) -> str:
    for line in text.splitlines():
        line = line.strip()
        if line.startswith(prefix):
            return line[len(prefix):]
    raise AssertionError(f"'{prefix}' not found in output:\n{text}")


# -------------------------------------------------------------
# Violations
# -------------------------------------------------------------


@pytest.mark.integration
def test_violations_agrees_with_cpp_voltage_check(tests_data_dir):
    """Agree with checkVoltageViolations(); its bool is the only oracle."""
    r = run_inline(
        """
        from gridpack import Session, PowerFlow
        with Session() as s:
            pf = PowerFlow(s, "input_14.xml", suppress_output=True)
            result = pf.solve()
            # Deliberately not bound to a local: a surviving reference to
            # the pybind11 app outlives the wrapper and aborts at MPI
            # finalize (MPI_Iprobe on MPI_COMM_NULL).
            for lo, hi in [(0.9, 1.1), (1.02, 1.06), (1.0, 1.2), (1.05, 1.1)]:
                pf._pfapp.setVoltageLimits(lo, hi)
                cpp_ok = pf._pfapp.checkVoltageViolations()  # True == none
                pf._pfapp.clearVoltageViolations()
                py = result.violations(min_voltage=lo, max_voltage=hi)["voltage"]
                print(f"AGREE={cpp_ok == (len(py) == 0)} "
                      f"band={lo},{hi} n={len(py)}")
        """,
        cwd=tests_data_dir,
    )
    assert r.returncode == 0, f"stdout:\n{r.stdout}\nstderr:\n{r.stderr}"
    agrees = [ln for ln in r.stdout.splitlines() if ln.strip().startswith("AGREE=")]
    assert len(agrees) == 4, f"expected 4 bands:\n{r.stdout}"
    assert all("AGREE=True" in ln for ln in agrees), f"disagreement:\n{r.stdout}"
    # A band that flags nothing would make the agreement vacuous.
    assert any("n=6" in ln for ln in agrees), f"no band found violations:\n{r.stdout}"


@pytest.mark.integration
def test_violations_reports_unrated_branches(tests_data_dir):
    """IEEE14.raw has no ratings, which must not read as "no overloads"."""
    r = run_inline(
        """
        from gridpack import Session, PowerFlow
        with Session() as s:
            pf = PowerFlow(s, "input_14.xml", suppress_output=True)
            v = pf.solve().violations()
            print(f"N_BRANCHES={v['n_branches']}")
            print(f"UNRATED={v['unrated_branches']}")
            print(f"OVERLOAD={len(v['overload'])}")
        """,
        cwd=tests_data_dir,
    )
    assert r.returncode == 0, f"stdout:\n{r.stdout}\nstderr:\n{r.stderr}"
    assert "N_BRANCHES=20" in r.stdout
    assert "UNRATED=20" in r.stdout
    assert "OVERLOAD=0" in r.stdout


@pytest.mark.integration
def test_violations_detects_overload(tests_data_dir, tmp_path, rated_raw):
    """Halved ratings put several branches over 100%; the C++ check agrees."""
    d = _rated_case(tests_data_dir, tmp_path, rated_raw, scale=0.5)
    r = run_inline(
        """
        from gridpack import Session, PowerFlow
        with Session() as s:
            pf = PowerFlow(s, "case.xml", suppress_output=True)
            result = pf.solve()
            v = result.violations()
            print(f"UNRATED={v['unrated_branches']}")
            print(f"OVERLOAD={len(v['overload'])}")
            print(f"WORST={max(o['loadingPercent'] for o in v['overload']):.1f}")
            print(f"AGREE={pf._pfapp.checkLineOverloadViolations() is False}")
        """,
        cwd=d,
    )
    assert r.returncode == 0, f"stdout:\n{r.stdout}\nstderr:\n{r.stderr}"
    assert "UNRATED=0" in r.stdout
    assert "AGREE=True" in r.stdout
    n = int(_extract(r.stdout, "OVERLOAD="))
    assert n >= 3, f"expected several overloads, got {n}:\n{r.stdout}"
    assert float(_extract(r.stdout, "WORST=")) > 150.0


@pytest.mark.integration
def test_violations_threshold_and_unrated_on_rated_case(
    tests_data_dir, tmp_path, rated_raw
):
    """At rating, nothing is overloaded; at 80% several branches are."""
    d = _rated_case(tests_data_dir, tmp_path, rated_raw)
    r = run_inline(
        """
        from gridpack import Session, PowerFlow
        with Session() as s:
            pf = PowerFlow(s, "case.xml", suppress_output=True)
            result = pf.solve()
            print(f"UNRATED={result.violations()['unrated_branches']}")
            print(f"AT100={len(result.violations()['overload'])}")
            v80 = result.violations(overload_threshold=80.0)
            print(f"AT80={len(v80['overload'])}")
        """,
        cwd=d,
    )
    assert r.returncode == 0, f"stdout:\n{r.stdout}\nstderr:\n{r.stderr}"
    assert "UNRATED=0" in r.stdout
    assert "AT100=0" in r.stdout
    assert int(_extract(r.stdout, "AT80=")) >= 3


# -------------------------------------------------------------
# Bulk export
# -------------------------------------------------------------


@pytest.mark.integration
def test_to_csv_writes_every_table(tests_data_dir, tmp_path):
    """Each table exports with its full header; the legacy view still works."""
    r = run_inline(
        f"""
        from gridpack import Session, PowerFlow
        out = {str(tmp_path)!r}
        with Session() as s:
            pf = PowerFlow(s, "input_14.xml", suppress_output=True)
            result = pf.solve()
            for t in ("buses", "branches", "generators"):
                result.to_csv(f"{{out}}/{{t}}.csv", table=t)
            result.to_csv(f"{{out}}/legacy.csv", [1, 2, 3])
            print("WROTE=ok")
        """,
        cwd=tests_data_dir,
    )
    assert r.returncode == 0, f"stdout:\n{r.stdout}\nstderr:\n{r.stderr}"
    assert "WROTE=ok" in r.stdout

    header = (tmp_path / "buses.csv").read_text().splitlines()
    assert header[0].startswith("busId,type,area,zone,baseKV,voltage,angle")
    assert len(header) == 15                       # 14 buses + header
    assert len((tmp_path / "branches.csv").read_text().splitlines()) == 21
    assert (tmp_path / "generators.csv").read_text().splitlines()[0].startswith(
        "busId,genId,pGen,qGen"
    )
    # The narrow per-bus view is unchanged.
    legacy = (tmp_path / "legacy.csv").read_text().splitlines()
    assert legacy[0] == "bus,vmag,vangle"
    assert len(legacy) == 4


@pytest.mark.integration
def test_to_records_rejects_unknown_table(tests_data_dir):
    r = run_inline(
        """
        from gridpack import Session, PowerFlow
        with Session() as s:
            pf = PowerFlow(s, "input_14.xml", suppress_output=True)
            result = pf.solve()
            try:
                result.to_records(table="nope")
            except ValueError as e:
                print(f"RAISED={'nope' in str(e)}")
        """,
        cwd=tests_data_dir,
    )
    assert r.returncode == 0, f"stdout:\n{r.stdout}\nstderr:\n{r.stderr}"
    assert "RAISED=True" in r.stdout
