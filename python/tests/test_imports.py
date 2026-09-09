# -------------------------------------------------------------
# file: python/tests/test_imports.py
# -------------------------------------------------------------
# Pure-Python import surface checks.  These run even when GridPACK
# itself isn't functional -- as long as the compiled extension can
# be loaded (which it must be for the package to install), the
# high-level API surface is present and typed correctly.
# -------------------------------------------------------------

import re
from pathlib import Path
from types import ModuleType

import pytest

from .conftest import run_inline


def test_top_level_symbols():
    from gridpack import (
        Session,
        PowerFlow,
        PowerFlowResult,
        DynamicSim,
        DynamicSimStepper,
        DSFResult,
        StateEstimation,
        StateEstimationResult,
        Hadrec,
    )
    for cls in (
        Session,
        PowerFlow,
        PowerFlowResult,
        DynamicSim,
        DynamicSimStepper,
        DSFResult,
        StateEstimation,
        StateEstimationResult,
        Hadrec,
    ):
        assert isinstance(cls, type)


def test_submodule_shims():
    import gridpack.hadrec
    import gridpack.dynamic_simulation
    import gridpack.powerflow
    import gridpack.state_estimation
    import gridpack.emt

    assert isinstance(gridpack.hadrec.Module, type)
    assert isinstance(gridpack.hadrec.Action, type)
    assert isinstance(gridpack.dynamic_simulation.Event, type)
    assert isinstance(gridpack.dynamic_simulation.EventVector, type)
    assert isinstance(gridpack.powerflow.Powerflow, type)
    assert isinstance(gridpack.state_estimation.SEApp, type)
    assert isinstance(gridpack.state_estimation.Measurement, type)


def test_hadrec_action_builders():
    from gridpack.hadrec import (
        make_load_shed_action,
        make_line_trip_action,
        make_generator_trip_action,
        ACT_LOAD_SHED,
        ACT_LINE_TRIP,
        ACT_GEN_TRIP,
    )

    a = make_load_shed_action(bus_number=5, percentage=-0.2)
    assert a.actiontype == ACT_LOAD_SHED
    assert a.bus_number == 5
    assert a.componentID == "1"
    assert a.percentage == -0.2

    b = make_line_trip_action(from_bus=6, to_bus=7)
    assert b.actiontype == ACT_LINE_TRIP
    assert b.brch_from_bus_number == 6
    assert b.brch_to_bus_number == 7

    c = make_generator_trip_action(bus_number=3, gen_id="1")
    assert c.actiontype == ACT_GEN_TRIP
    assert c.bus_number == 3


# pyproject is the only place the release number is written by hand; the
# metadata and the CLI both derive from it, and these pin them to it.
_PYPROJECT = Path(__file__).resolve().parents[1] / "pyproject.toml"


def _pyproject_version():
    """Our own pyproject version, or None if that is not what sits there.

    The tests ship without the source tree: in the Docker image the directory
    above them holds the old setup.py-era pyproject, which has no [project].
    """
    try:
        txt = _PYPROJECT.read_text()
    except OSError:
        return None
    if "[project]" not in txt:
        return None
    m = re.search(r'^version = "([^"]+)"', txt, re.M)
    assert m, "no version in %s" % _PYPROJECT
    return m.group(1)


def test_metadata_matches_pyproject():
    from importlib.metadata import version
    expected = _pyproject_version()
    if expected is None:
        pytest.skip("%s is not this package's pyproject" % _PYPROJECT)
    assert version("gridpack") == expected


def test_cli_version_comes_from_package(capsys):
    """A separate CLI version number could not be kept in step, so there is none."""
    from importlib.metadata import version
    from gridpack.cli.main import main
    with pytest.raises(SystemExit) as excinfo:
        main(["--version"])
    assert excinfo.value.code == 0
    assert capsys.readouterr().out.strip() == "gridpack %s" % version("gridpack")


# -------------------------------------------------------------
# Curated namespace
# -------------------------------------------------------------
# `from ._gridpack import *` used to pre-bind the pybind11 submodules as
# attributes of this package, so CPython skipped importing the same-named
# Python shims and gridpack.dynamic_simulation.DynamicSim raised
# AttributeError.  The shim checks run in a subprocess because
# test_submodule_shims above imports those shims by statement, which
# populates sys.modules and hides the fault entirely.

_SHIM_DRIVER = """
    import sys
    import gridpack                       # bare import only

    for name in gridpack.COMPAT_MODULES:
        attr = getattr(gridpack, name)
        assert attr is sys.modules.get("gridpack." + name), name
        assert attr.__file__.endswith(name + ".py"), (name, attr.__file__)

    # dynamic_simulation.py re-exports these so the compat name reaches the
    # high-level classes; emt.py is where the EMT wrapper will land.
    assert gridpack.dynamic_simulation.DynamicSim is gridpack.DynamicSim
    print("SHIMS OK")
"""


def test_shims_resolve_to_python_modules_after_bare_import(tmp_path):
    r = run_inline(_SHIM_DRIVER, cwd=tmp_path)
    assert r.returncode == 0, f"stdout:\n{r.stdout}\nstderr:\n{r.stderr}"
    assert "SHIMS OK" in r.stdout


def test_all_names_are_importable():
    import gridpack
    assert not [n for n in gridpack.__all__ if not hasattr(gridpack, n)]


def test_all_keeps_the_low_level_types():
    """The scripts under python/src drive GridPACK through these."""
    import gridpack
    for n in ("CoarseTimer", "Communicator", "Configuration",
              "ConfigurationCursor", "Environment", "NoPrint",
              "TaskCounter", "TaskManager"):
        assert n in gridpack.__all__, n
        assert isinstance(getattr(gridpack, n), type), n


def test_dir_hides_implementation_modules():
    import gridpack
    d = dir(gridpack)
    for n in ("session", "results", "exceptions", "dynamic_sim", "contingency"):
        assert n not in d, n
    for n in gridpack.COMPAT_MODULES:
        assert n in d, n


def test_import_star_binds_no_modules():
    ns = {}
    exec("from gridpack import *", ns)
    assert not [k for k, v in ns.items() if isinstance(v, ModuleType)]


def test_unknown_attribute_names_the_public_api():
    import gridpack
    with pytest.raises(AttributeError, match="the public API is"):
        gridpack.NoSuchThing


# -------------------------------------------------------------
# EMT is experimental
# -------------------------------------------------------------

def test_emt_use_warns_experimental():
    import gridpack
    with pytest.warns(FutureWarning, match="experimental"):
        gridpack.emt.EMT


_QUIET_DRIVER = """
    import warnings
    warnings.simplefilter("error", FutureWarning)
    import gridpack               # emt is imported eagerly; must stay quiet
    print("QUIET OK")
"""


def test_plain_import_does_not_warn(tmp_path):
    """The notice has to fire on use, or it is noise on every import."""
    r = run_inline(_QUIET_DRIVER, cwd=tmp_path)
    assert r.returncode == 0, f"stdout:\n{r.stdout}\nstderr:\n{r.stderr}"
    assert "QUIET OK" in r.stdout
