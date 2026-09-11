# -------------------------------------------------------------
# file: python/tests/test_contingency.py
# -------------------------------------------------------------
# Tests for gridpack.ContingencyAnalysis and the contingency-list parser.
# -------------------------------------------------------------

import re

import pytest

from .conftest import CA_EXPECTED_STATUS, _data_sets, run_inline

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


def test_parse_tolerates_indented_declaration(tmp_path):
    """ca.x reads the shipped 118/polish/euro lists; ElementTree would not.

    Those three indent the <?xml?> declaration, and a declaration that is
    not at byte 0 is a hard parse error.
    """
    p = tmp_path / "ctg.xml"
    p.write_text(
        '  \n<?xml version="1.0" encoding="utf-8"?>\n'
        "<ContingencyList><Contingency_analysis><Contingencies>\n"
        """      <Contingency>
        <contingencyType>Line</contingencyType>
        <contingencyName>N1</contingencyName>
        <contingencyLineBuses>2 3</contingencyLineBuses>
        <contingencyLineNames>BL</contingencyLineNames>
      </Contingency>"""
        "\n</Contingencies></Contingency_analysis></ContingencyList>\n"
    )
    (c,) = parse_contingency_list(p)
    assert c.name == "N1" and c.from_buses == [2]


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
    # Neither is a solver failure, so neither reads as DIVERGENT.
    (dict(found=True, converged=False, islanded=True), "ISLANDED"),
    (dict(found=True, converged=False, slack_overload=True), "SLACK OVERLOAD"),
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


def test_result_reports_slack_overload_like_ca_x():
    """ca.x wording: the solve succeeded but the slack exceeded Pmax."""
    r = ContingencyResult("C", found=True, converged=False, slack_overload=True)
    assert r.report_lines == [
        "Insufficient generation capacity for contingency C"]
    assert not r.ok


def test_result_reports_islanding():
    r = ContingencyResult("C", found=True, converged=False, islanded=True)
    assert r.report_lines == ["Islanding detected for contingency C"]
    assert not r.ok


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
        results = ca.gather()
        # Rank 0 only: both ranks writing the same lines into one pipe
        # interleaved mid-line, which _statuses then mis-parsed.
        if s.rank == 0:
            for r in results:
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
    # One rank prints; two ranks sharing the pipe used to interleave mid-line.
    assert r.stdout.count("RESULT ") == len(CA_EXPECTED_STATUS), r.stdout
    assert _statuses(r.stdout) == CA_EXPECTED_STATUS
    # gather() sorts by name, so both ranks report the same order regardless
    # of which rank drew which task.
    assert list(_statuses(r.stdout)) == sorted(CA_EXPECTED_STATUS)


@pytest.mark.integration
def test_contingency_analysis_warns_on_group_size(ca_case):
    """groupSize > 1 is silently downgraded, so it has to say so."""
    cfg = ca_case / "input_14.xml"
    # Strip then insert: the shared XML has carried groupSize and then not.
    src = re.sub(r"\s*<groupSize>\d+</groupSize>", "", cfg.read_text())
    src = src.replace("<Contingency_analysis>",
                      "<Contingency_analysis>\n    <groupSize>2</groupSize>", 1)
    assert "<groupSize>2</groupSize>" in src, "config edit did nothing"
    cfg.write_text(src)
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
    # 0.98/1.02, not tighter: ignore_voltage_violations() exempts whatever
    # the base case already violates, so a very tight band reports fewer
    # bus violations, not more.
    assert _statuses(r.stdout)["LINE_2_3"] == "BUS+BRANCH VIOLATION"


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
        "Branch Violation for contingency LINE_2_3")
    assert (ca_case / "GEN_2.out").read_text().rstrip().endswith(
        "Divergent for contingency GEN_2")


@pytest.mark.integration
def test_slack_overload_is_not_reported_as_converged(ca_case):
    """A solve that only balances by overdrawing the slack is not a solution.

    IEEE14.raw has less generation headroom than the IEEE14_ca.raw the case
    normally uses, and ca.x reports "Insufficient generation capacity" for
    all three line outages on it.  Without checkSlackCapacity they came
    back OK, which is the worst possible answer: a masked failure.
    """
    import shutil
    shutil.copy(_data_sets() / "raw/IEEE14.raw", ca_case / "IEEE14_ca.raw")
    r = run_inline(_DRIVER, cwd=ca_case)
    assert r.returncode == 0, r.stderr[-2000:]
    got = _statuses(r.stdout)
    assert [got[n] for n in ("LINE_2_3", "LINE_6_13", "LINE_13_14")] == \
        ["SLACK OVERLOAD"] * 3, got


@pytest.mark.integration
def test_118_status_split_matches_ca_x(tmp_path):
    """Pin the status split on IEEE-118 to what ca.x reports on the same run.

    The 14-bus case splits into no islands, so this is the only coverage of
    the pre-solve island check.  Counts come from ca.x's own
    ca_results_convergence.csv: 160 solved, 17 SLACK_OVERLOAD, 2 ISLANDED.
    """
    import shutil
    data = _data_sets()
    src = [data / "input/ca/input_118.xml",
           data / "contingencies/contingencies_118.xml",
           data / "raw/IEEE118.raw"]
    for f in src:
        if not f.exists():
            pytest.skip("missing %s" % f)
        shutil.copy(f, tmp_path / f.name)

    r = run_inline("""
        from collections import Counter
        from gridpack import Session, ContingencyAnalysis
        with Session() as s:
            ca = ContingencyAnalysis(s, "input_118.xml",
                                     print_calc_files=False,
                                     suppress_output=True)
            ca.run()
            c = Counter(x.status for x in ca.gather())
            if s.rank == 0:
                for k in sorted(c):
                    print("COUNT %s=%d" % (k, c[k]))
    """, cwd=tmp_path, timeout=600)
    assert r.returncode == 0, r.stderr[-2000:]
    counts = dict(
        (kv.split("=")[0], int(kv.split("=")[1]))
        for kv in (l.split(" ", 1)[1] for l in r.stdout.splitlines()
                   if l.startswith("COUNT ")))
    assert counts.get("ISLANDED") == 2, counts
    assert counts.get("SLACK OVERLOAD") == 17, counts
    # The rest solved; ca.x finds 5 bus and 160 branch violations among them.
    solved = sum(v for k, v in counts.items()
                 if k not in ("ISLANDED", "SLACK OVERLOAD"))
    assert solved == 160, counts
    assert counts.get("BUS+BRANCH VIOLATION") == 5, counts


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
