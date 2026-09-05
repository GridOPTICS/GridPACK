# -------------------------------------------------------------
# file: python/tests/test_contingency.py
# -------------------------------------------------------------
# Tests for gridpack.ContingencyAnalysis and the contingency-list parser.
# -------------------------------------------------------------

import pytest

from .conftest import CA_EXPECTED_STATUS, run_inline

from gridpack.contingency import (
    Contingency,
    ContingencyResult,
    parse_contingency_list,
)


# -------------------------------------------------------------
# Contingency-list parsing
# -------------------------------------------------------------

def _write_list(tmp_path, body: str):
    p = tmp_path / "ctg.xml"
    p.write_text(
        '<?xml version="1.0" encoding="utf-8"?>\n'
        "<ContingencyList><Contingency_analysis><Contingencies>\n"
        + body +
        "\n</Contingencies></Contingency_analysis></ContingencyList>\n"
    )
    return p


def test_parse_line_contingency_pairs_buses(tmp_path):
    """Line buses are a flat from/to list, one pair per circuit id."""
    p = _write_list(tmp_path, """
      <Contingency>
        <contingencyType>Line</contingencyType>
        <contingencyName>N2</contingencyName>
        <contingencyLineBuses>2 3 6 13</contingencyLineBuses>
        <contingencyLineNames>BL B2</contingencyLineNames>
      </Contingency>""")
    (c,) = parse_contingency_list(p)
    assert c.kind == "Line" and c.is_line
    assert c.name == "N2"
    assert c.from_buses == [2, 6]
    assert c.to_buses == [3, 13]
    assert c.circuit_ids == ["BL", "B2"]


def test_parse_generator_contingency(tmp_path):
    p = _write_list(tmp_path, """
      <Contingency>
        <contingencyType>Generator</contingencyType>
        <contingencyName>GEN_2</contingencyName>
        <contingencyBuses>2 3</contingencyBuses>
        <contingencyGenerators>1 2</contingencyGenerators>
      </Contingency>""")
    (c,) = parse_contingency_list(p)
    assert not c.is_line
    assert c.buses == [2, 3]
    # PSS/E ids are 2 chars wide; "1" must become "1 ", not "1".
    assert c.generator_ids == ["1 ", "2 "]


def test_parse_pads_and_truncates_circuit_ids(tmp_path):
    """The C++ side compares fixed-width ids, so padding is not cosmetic."""
    p = _write_list(tmp_path, """
      <Contingency>
        <contingencyType>Line</contingencyType>
        <contingencyName>W</contingencyName>
        <contingencyLineBuses>1 2 2 3</contingencyLineBuses>
        <contingencyLineNames>1 LONG</contingencyLineNames>
      </Contingency>""")
    (c,) = parse_contingency_list(p)
    assert c.circuit_ids == ["1 ", "LO"]


def test_parse_rejects_unknown_type(tmp_path):
    """A typo in contingencyType silently dropped the contingency before."""
    p = _write_list(tmp_path, """
      <Contingency>
        <contingencyType>Transformer</contingencyType>
        <contingencyName>X</contingencyName>
      </Contingency>""")
    with pytest.raises(ValueError, match="unknown contingencyType"):
        parse_contingency_list(p)


def test_parse_rejects_short_bus_list(tmp_path):
    """Two names but one bus pair used to IndexError deep in the loop."""
    p = _write_list(tmp_path, """
      <Contingency>
        <contingencyType>Line</contingencyType>
        <contingencyName>X</contingencyName>
        <contingencyLineBuses>1 2</contingencyLineBuses>
        <contingencyLineNames>BL B2</contingencyLineNames>
      </Contingency>""")
    with pytest.raises(ValueError, match="2 buses per name"):
        parse_contingency_list(p)


def test_parse_empty_list(tmp_path):
    assert parse_contingency_list(_write_list(tmp_path, "")) == []


# -------------------------------------------------------------
# ContingencyResult reporting
# -------------------------------------------------------------

@pytest.mark.parametrize("kwargs,status", [
    (dict(found=True, converged=True), "OK"),
    (dict(found=True, converged=True, voltage_ok=False), "BUS VIOLATION"),
    (dict(found=True, converged=True, overload_ok=False), "BRANCH VIOLATION"),
    (dict(found=True, converged=True, voltage_ok=False, overload_ok=False),
     "BUS+BRANCH VIOLATION"),
    (dict(found=True, converged=False), "DIVERGENT"),
    # Not-found is a config error, not a network finding -- keep it distinct.
    (dict(found=False, converged=False), "NOT FOUND"),
])
def test_result_status(kwargs, status):
    assert ContingencyResult(name="C", **kwargs).status == status


def test_result_ok_only_when_clean():
    assert ContingencyResult("C", found=True, converged=True).ok
    assert not ContingencyResult("C", True, True, voltage_ok=False).ok
    assert not ContingencyResult("C", found=False, converged=False).ok


def test_result_reports_both_violations():
    """A doubly-violating case gets one line per violation, as ca.x does."""
    r = ContingencyResult("C", True, True, voltage_ok=False, overload_ok=False)
    assert r.report_lines == ["Bus Violation for contingency C",
                             "Branch Violation for contingency C"]


def test_result_reports_no_violation():
    r = ContingencyResult("C", found=True, converged=True)
    assert r.report_lines == ["No violation for contingency C"]


def test_line_contingency_to_pybind_sets_branch_type():
    c = Contingency(name="L", kind="Line", from_buses=[2], to_buses=[3],
                    circuit_ids=["BL"])
    ctg = c.to_pybind()
    assert ctg.p_name == "L"
    assert ctg.p_type == 1
    assert list(ctg.p_ckt) == ["BL"]
    assert list(ctg.p_saveLineStatus) == [1]


def test_generator_contingency_to_pybind_sets_generator_type():
    c = Contingency(name="G", kind="Generator", buses=[2], generator_ids=["1 "])
    ctg = c.to_pybind()
    assert ctg.p_type == 0
    assert list(ctg.p_busid) == [2]
    assert list(ctg.p_saveGenStatus) == [1]


# -------------------------------------------------------------
# Live runs
# -------------------------------------------------------------

_DRIVER = """
    from gridpack import Session, ContingencyAnalysis
    with Session() as s:
        ca = ContingencyAnalysis(s, "input_14.xml", print_calc_files=False,
                                 suppress_output=True)
        ca.run()
        for r in ca.gather():
            print("RESULT %s %s" % (r.name, r.status))
"""


def _statuses(stdout):
    out = {}
    for line in stdout.splitlines():
        if line.strip().startswith("RESULT "):
            _, name, status = line.strip().split(" ", 2)
            out[name] = status
    return out


@pytest.mark.integration
def test_contingency_analysis_serial(ca_case):
    r = run_inline(_DRIVER, cwd=ca_case)
    assert r.returncode == 0, f"stdout:\n{r.stdout}\nstderr:\n{r.stderr}"
    assert _statuses(r.stdout) == CA_EXPECTED_STATUS


@pytest.mark.integration
@pytest.mark.mpi
def test_contingency_analysis_under_mpi_matches_serial(ca_case,
                                                       require_mpiexec):
    """Each task needs its own network copy or the collective solve deadlocks.

    On the world communicator, ranks entered solve() on different
    contingencies and rank 1 reported "No reference bus found".
    """
    r = run_inline(_DRIVER, cwd=ca_case, mpi_np=2, timeout=180)
    assert r.returncode == 0, f"stdout:\n{r.stdout}\nstderr:\n{r.stderr}"
    assert "No reference bus found" not in r.stdout + r.stderr
    assert _statuses(r.stdout) == CA_EXPECTED_STATUS
    # gather() sorts by name, so both ranks report the same order regardless
    # of which rank drew which task.
    assert list(_statuses(r.stdout)) == sorted(CA_EXPECTED_STATUS)


@pytest.mark.integration
def test_contingency_analysis_warns_on_group_size(ca_case):
    """groupSize > 1 is silently downgraded, so it has to say so."""
    cfg = ca_case / "input_14.xml"
    cfg.write_text(cfg.read_text().replace("<groupSize>1</groupSize>",
                                           "<groupSize>2</groupSize>"))
    r = run_inline(_DRIVER, cwd=ca_case)
    assert r.returncode == 0, r.stderr[-2000:]
    assert "groupSize=2 is ignored" in r.stderr
    assert _statuses(r.stdout) == CA_EXPECTED_STATUS


@pytest.mark.integration
def test_contingency_analysis_honors_voltage_limits(ca_case):
    """A tight band must actually reach checkVoltageViolations."""
    r = run_inline("""
        from gridpack import Session, ContingencyAnalysis
        with Session() as s:
            ca = ContingencyAnalysis(s, "input_14.xml",
                                     voltage_limits=(0.98, 1.02),
                                     print_calc_files=False,
                                     suppress_output=True)
            ca.run()
            for r in ca.gather():
                print("RESULT %s %s" % (r.name, r.status))
    """, cwd=ca_case)
    assert r.returncode == 0, r.stderr[-2000:]
    assert _statuses(r.stdout)["LINE_2_3"] == "BUS VIOLATION"


@pytest.mark.integration
def test_contingency_analysis_writes_calc_files(ca_case):
    """suppress_output is deliberately absent -- it gags open()/print() too."""
    r = run_inline("""
        from gridpack import Session, ContingencyAnalysis
        with Session() as s:
            ContingencyAnalysis(s, "input_14.xml").run()
    """, cwd=ca_case)
    assert r.returncode == 0, r.stderr[-2000:]
    assert (ca_case / "LINE_2_3.out").read_text().rstrip().endswith(
        "No violation for contingency LINE_2_3")
    assert (ca_case / "GEN_2.out").read_text().rstrip().endswith(
        "Divergent for contingency GEN_2")


@pytest.mark.integration
def test_contingency_analysis_closes_with_session(ca_case):
    """The wrapper must be registered, or PFAppModule outlives MPI."""
    r = run_inline("""
        from gridpack import Session, ContingencyAnalysis
        s = Session()
        ca = ContingencyAnalysis(s, "input_14.xml", print_calc_files=False,
                                 suppress_output=True)
        s.close()
        print("CLOSED", ca._closed)
        try:
            ca.run()
        except RuntimeError as e:
            print("RAISED", e)
    """, cwd=ca_case)
    assert r.returncode == 0, f"stdout:\n{r.stdout}\nstderr:\n{r.stderr}"
    assert "CLOSED True" in r.stdout
    assert "RAISED ContingencyAnalysis is closed" in r.stdout


@pytest.mark.integration
def test_suppress_output_also_suppresses_calc_files(ca_case):
    """suppressOutput gags the C++ stream, so open() writes nothing.

    Pinned because printCalcFiles=true suggests files should appear.
    """
    r = run_inline("""
        from gridpack import Session, ContingencyAnalysis
        with Session() as s:
            ContingencyAnalysis(s, "input_14.xml",
                                suppress_output=True).run()
    """, cwd=ca_case)
    assert r.returncode == 0, r.stderr[-2000:]
    assert not list(ca_case.glob("*.out"))
