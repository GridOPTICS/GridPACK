# -------------------------------------------------------------
# file: python/tests/conftest.py
# -------------------------------------------------------------
# Shared fixtures and helpers for the gridpack pytest suite.
#
# Integration tests exercise the compiled pybind11 extension and are
# marked ``@pytest.mark.integration``.  Because :class:`gridpack.Session`
# is a process singleton (MPI/PETSc/GA can only be initialized once per
# process), integration tests that need a session run the payload in
# an *inner* Python subprocess and assert on stdout / exit code.  This
# also lets us cover the step-by-step DS path that terminates with
# ``os._exit(0)`` (upstream ``DSFullApp::~`` SEGV workaround).
# -------------------------------------------------------------

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import textwrap
from pathlib import Path
from typing import List, Optional


def pytest_configure(config):
    config.addinivalue_line(
        "markers",
        "integration: exercises the compiled gridpack._gridpack extension "
        "(needs GRIDPACK_DIR / a build)",
    )
    config.addinivalue_line(
        "markers",
        "mpi: needs mpiexec available on PATH",
    )


# -------------------------------------------------------------
# Helpers
# -------------------------------------------------------------


def _tests_data_dir() -> Path:
    """Directory that holds a copy of the standard 14-bus test data."""
    return Path("/tmp/pygridpack-test").resolve()


def _dsf_build_dir() -> Path:
    """Location of ``input_9b3g.xml`` (in the DSF build tree)."""
    return Path(
        "/Users/yousu.chen/Software/GridPACK/src/build/applications/"
        "dynamic_simulation_full_y"
    ).resolve()


def _se_data_dir() -> Path:
    return _tests_data_dir() / "se"


def run_inline(
    body: str,
    *,
    cwd: Optional[Path] = None,
    mpi_np: int = 0,
    timeout: int = 120,
) -> subprocess.CompletedProcess:
    """Execute ``body`` as an inner Python script and return the result.

    Parameters
    ----------
    body : str
        Python source to write to a temp file and run.
    cwd : Path, optional
        Working directory for the subprocess.
    mpi_np : int, optional
        If > 0, launch under ``mpiexec -np <mpi_np>``.
    timeout : int, optional
        Kill the subprocess after this many seconds.
    """
    import tempfile

    with tempfile.NamedTemporaryFile(
        "w", suffix=".py", delete=False
    ) as fh:
        fh.write(textwrap.dedent(body))
        script = fh.name
    try:
        cmd: List[str] = []
        if mpi_np > 0:
            cmd.extend(["mpiexec", "-np", str(mpi_np)])
        cmd.extend([sys.executable, "-u", script])
        return subprocess.run(
            cmd,
            cwd=str(cwd) if cwd else None,
            capture_output=True,
            text=True,
            timeout=timeout,
            check=False,
        )
    finally:
        try:
            os.unlink(script)
        except OSError:
            pass


# -------------------------------------------------------------
# Fixtures
# -------------------------------------------------------------


import pytest


@pytest.fixture(scope="session")
def tests_data_dir() -> Path:
    """Directory holding IEEE14 test inputs.  Skips the suite if absent."""
    d = _tests_data_dir()
    if not (d / "input_14.xml").exists() or not (d / "IEEE14.raw").exists():
        pytest.skip(
            f"missing test data under {d} (input_14.xml / IEEE14.raw). "
            "Copy them into place to run integration tests."
        )
    return d


@pytest.fixture(scope="session")
def rated_raw() -> Path:
    """The in-repo IEEE-14 with synthetic line ratings.

    IEEE14.raw has rateA <= 0 on every branch, so it cannot exercise
    overloads at all.
    """
    d = Path(__file__).resolve().parents[2] / "src/applications/data_sets/raw"
    f = d / "IEEE14_PTIv33_rated.raw"
    if not f.exists():
        pytest.skip(f"missing {f}")
    return f


@pytest.fixture(scope="session")
def se_data_dir() -> Path:
    """Directory holding the state-estimation IEEE14 inputs."""
    d = _se_data_dir()
    for name in ("input.xml", "IEEE14.raw", "IEEE14_meas.xml"):
        if not (d / name).exists():
            pytest.skip(
                f"missing SE test data at {d}/{name}. Copy from "
                "src/applications/modules/state_estimation/test/"
            )
    return d


# The shipped contingencies_14.xml names circuits B1/B2/B3, but every branch
# in IEEE14.raw is "BL", so every contingency there comes back not-found.
# These hit the real network: three branches plus one diverging generator.
_CA_CONTINGENCIES = """<?xml version="1.0" encoding="utf-8"?>
<ContingencyList>
  <Contingency_analysis>
    <Contingencies>
      <Contingency>
        <contingencyType>Line</contingencyType>
        <contingencyName>LINE_2_3</contingencyName>
        <contingencyLineBuses>2 3</contingencyLineBuses>
        <contingencyLineNames>BL</contingencyLineNames>
      </Contingency>
      <Contingency>
        <contingencyType>Line</contingencyType>
        <contingencyName>LINE_6_13</contingencyName>
        <contingencyLineBuses>6 13</contingencyLineBuses>
        <contingencyLineNames>BL</contingencyLineNames>
      </Contingency>
      <Contingency>
        <contingencyType>Line</contingencyType>
        <contingencyName>LINE_13_14</contingencyName>
        <contingencyLineBuses>13 14</contingencyLineBuses>
        <contingencyLineNames>BL</contingencyLineNames>
      </Contingency>
      <Contingency>
        <contingencyType>Generator</contingencyType>
        <contingencyName>GEN_2</contingencyName>
        <contingencyBuses>2</contingencyBuses>
        <contingencyGenerators>1</contingencyGenerators>
      </Contingency>
    </Contingencies>
  </Contingency_analysis>
</ContingencyList>
"""

# Statuses the four contingencies above produce at the default 0.9-1.1 band.
CA_EXPECTED_STATUS = {
    "LINE_2_3": "OK",
    "LINE_6_13": "OK",
    "LINE_13_14": "OK",
    "GEN_2": "DIVERGENT",
}


@pytest.fixture
def ca_case(tmp_path) -> Path:
    """A runnable CA case directory, returned as its path.

    Switches the in-repo input off FullBranchN1 (unimplemented) onto a list.
    """
    data = Path(__file__).resolve().parents[2] / "src/applications/data_sets"
    cfg_src = data / "input/ca/input_14.xml"
    raw_src = data / "raw/IEEE14.raw"
    for f in (cfg_src, raw_src):
        if not f.exists():
            pytest.skip(f"missing {f}")

    cfg = cfg_src.read_text()
    cfg = cfg.replace("<FullBranchN1>true</FullBranchN1>",
                      "<FullBranchN1>false</FullBranchN1>")
    cfg = cfg.replace("<FullGeneratorN1>true</FullGeneratorN1>",
                      "<FullGeneratorN1>false</FullGeneratorN1>")
    marker = "<!-- Option 2: Auto-generate N-1 contingencies from network -->"
    assert marker in cfg, "CA input format changed; contingencyList not wired"
    cfg = cfg.replace(
        marker, "<contingencyList>contingencies.xml</contingencyList>")

    (tmp_path / "input_14.xml").write_text(cfg)
    (tmp_path / "contingencies.xml").write_text(_CA_CONTINGENCIES)
    shutil.copy(raw_src, tmp_path / "IEEE14_ca.raw")
    return tmp_path


@pytest.fixture(scope="session")
def dsf_build_dir() -> Path:
    """Directory holding ``input_9b3g.xml`` (in the DSF build tree)."""
    d = _dsf_build_dir()
    if not (d / "input_9b3g.xml").exists():
        pytest.skip(
            f"missing {d}/input_9b3g.xml; build the DSF example first."
        )
    return d


@pytest.fixture(scope="session")
def has_mpiexec() -> bool:
    return shutil.which("mpiexec") is not None


@pytest.fixture
def require_mpiexec(has_mpiexec) -> None:
    if not has_mpiexec:
        import pytest as _pytest
        _pytest.skip("mpiexec not on PATH")
