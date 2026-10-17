# -------------------------------------------------------------
# file: python/tests/test_case.py
# -------------------------------------------------------------
# Case: load, edit and stage an input XML.
# -------------------------------------------------------------

import xml.etree.ElementTree as ET

import pytest

from gridpack import Case
from .conftest import _data_sets, run_inline


def _written(xml):
    return ET.parse(xml).getroot()


def test_reads_what_boost_accepts(tmp_path):
    # <!-ksp_view> is not a comment; Boost skips it, ElementTree rejects it.
    case = Case.from_data_sets("input/ds/input_9b3g.xml")
    opts = case.get("Powerflow.LinearSolver.PETScOptions")
    assert "-ksp_type richardson" in opts and "ksp_view" not in opts
    root = _written(case.write(tmp_path))
    assert root.find("Powerflow/LinearSolver/PETScOptions").text.split() \
        == opts.split()


def test_write_copies_references_as_bare_names(tmp_path):
    case = Case.from_data_sets("input/ds/input_9b3g.xml")
    assert case.files() == {
        "Powerflow.networkConfiguration": str(_data_sets() / "raw/9b3g.raw"),
        "Dynamic_simulation.generatorParameters":
            str(_data_sets() / "dyr/9b3g.dyr"),
    }
    xml = case.write(tmp_path)
    assert sorted(p.name for p in tmp_path.iterdir()) \
        == ["9b3g.dyr", "9b3g.raw", "input_9b3g.xml"]
    root = _written(xml)
    assert root.find("Powerflow/networkConfiguration").text == "9b3g.raw"
    # Outputs are not references.
    assert root.find("Dynamic_simulation/generatorWatchFileName").text.strip() \
        == "gen_watch_9b3g.csv"


def test_set_creates_missing_elements_and_formats_values():
    case = Case.from_data_sets("input/ca/input_14.xml")
    case.set("Powerflow.tolerance", 1e-8)
    case.set("Contingency_analysis.FullBranchN1", False)
    case.set("Configuration.New.Deeper.leaf", 3)
    assert case.get("Powerflow.tolerance") == "1e-08"
    assert case.get("Contingency_analysis.FullBranchN1") == "false"
    assert case.get("New.Deeper.leaf") == "3"
    assert case.get("New.missing", "x") == "x"
    with pytest.raises(ValueError):
        case.set("Powerflow", 1)


def test_set_resolves_a_reference_in_data_sets():
    case = Case.from_data_sets("input/powerflow/input_14.xml")
    case.set("Powerflow.networkConfiguration", "IEEE14_PTIv33_rated.raw")
    assert case.get("Powerflow.networkConfiguration") \
        == str(_data_sets() / "raw/IEEE14_PTIv33_rated.raw")


def test_missing_reference_fails_at_write_unless_deferred(tmp_path):
    case = Case.from_data_sets("input/ca/input_14.xml")
    case.set("Contingency_analysis.contingencyList", "later.xml")
    with pytest.raises(FileNotFoundError, match="contingencyList"):
        case.write(tmp_path)
    root = _written(case.write(tmp_path, deferred=["later.xml"]))
    assert root.find("Contingency_analysis/contingencyList").text == "later.xml"


def test_include_keeps_the_other_cases_references(tmp_path):
    case = Case.from_data_sets("input/powerflow/input_14.xml")
    case.include("input/se/input_14.xml", "State_estimation")
    case.write(tmp_path)
    assert (tmp_path / "IEEE14_meas.xml").is_file()
    case.include("input/se/input_14.xml", "State_estimation")
    assert len(_written(case.write(tmp_path)).findall("State_estimation")) == 1
    with pytest.raises(KeyError):
        case.include("input/se/input_14.xml", "Dynamic_simulation")


def test_set_xml_replaces_content():
    case = Case.from_data_sets("input/ds/input_9b3g.xml")
    case.set_xml("Dynamic_simulation.observations",
                 "<observation><type>bus</type><busID>1</busID></observation>")
    assert case.get("Dynamic_simulation.observations.observation.busID") == "1"


def test_name_clash_and_bad_xml_are_reported(tmp_path):
    case = Case.from_data_sets("input/powerflow/input_14.xml")
    (tmp_path / "IEEE14.raw").write_text("")
    case.set("Dynamic_simulation.generatorParameters", tmp_path / "IEEE14.raw")
    with pytest.raises(ValueError, match="two different files named IEEE14.raw"):
        case.write(tmp_path / "out")

    broken = tmp_path / "broken.xml"
    broken.write_text("<Configuration><a>1</b></Configuration>")
    with pytest.raises(ValueError, match="broken.xml"):
        Case(broken)


@pytest.mark.integration
def test_edits_reach_the_solver(tmp_path):
    case = Case.from_data_sets("input/powerflow/input_14.xml")
    case.set("Powerflow.maxIteration", 1)
    case.set("Powerflow.tolerance", 1e-12)
    case.write(tmp_path)
    r = run_inline(
        """
        from gridpack import Session, PowerFlow
        with Session() as s:
            res = PowerFlow(s, "input_14.xml", suppress_output=True).solve(
                strict=False)
            print("CONVERGED=%s ITER=%s" % (res.converged, res.iterations))
        """,
        cwd=tmp_path,
    )
    assert r.returncode == 0, f"stdout:\n{r.stdout}\nstderr:\n{r.stderr}"
    assert "CONVERGED=False ITER=1" in r.stdout
