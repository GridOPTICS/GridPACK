# -------------------------------------------------------------
# file: python/tests/test_cli.py
# -------------------------------------------------------------
# CLI smoke tests: --help, --version.  A live PF run through the CLI
# is exercised by test_powerflow's integration case.
# -------------------------------------------------------------

import shutil
import subprocess
from pathlib import Path

import pytest


def _has_gridpack_cli():
    return shutil.which("gridpack") is not None


@pytest.mark.skipif(not _has_gridpack_cli(),
                    reason="gridpack console script not installed")
def test_cli_help():
    r = subprocess.run(
        ["gridpack", "--help"], capture_output=True, text=True, check=False,
    )
    assert r.returncode == 0
    for cmd in ("powerflow", "dsf", "se", "hadrec", "emt", "ca"):
        assert cmd in r.stdout, f"'{cmd}' subcommand missing from --help"


@pytest.mark.skipif(not _has_gridpack_cli(),
                    reason="gridpack console script not installed")
def test_cli_version():
    r = subprocess.run(
        ["gridpack", "--version"], capture_output=True, text=True, check=False,
    )
    # Some argparse configurations write --version to stdout, some to
    # stderr; accept either.
    out = (r.stdout or "") + (r.stderr or "")
    assert r.returncode == 0
    assert out.strip(), "gridpack --version produced no output"


# -------------------------------------------------------------
# Exit codes
# -------------------------------------------------------------
# 0 success / 1 unexpected error / 2 bad usage or config / 3 diverged.
# Scripts branch on these, so each one is pinned here.

_CLI = pytest.mark.skipif(not _has_gridpack_cli(),
                          reason="gridpack console script not installed")

_DATA = Path(__file__).resolve().parents[2] / "src/applications/data_sets"
_CA_INPUT = _DATA / "input/ca/input_14.xml"
_RAW_14 = _DATA / "raw/IEEE14.raw"


def _run(*argv, cwd=None, np=0, env=None):
    cmd = (["mpiexec", "-np", str(np)] if np else []) + ["gridpack", *argv]
    return subprocess.run(cmd, cwd=str(cwd) if cwd else None,
                          capture_output=True, text=True, check=False,
                          timeout=300, env=env)


def _capped(tests_data_dir, tmp_path, iterations=2):
    """Copy the 14-bus case with maxIteration too low to converge."""
    src = (tests_data_dir / "input_14.xml").read_text()
    dst = tmp_path / "capped.xml"
    dst.write_text(src.replace("<maxIteration>50</maxIteration>",
                               f"<maxIteration>{iterations}</maxIteration>"))
    shutil.copy(tests_data_dir / "IEEE14.raw", tmp_path / "IEEE14.raw")
    return dst


@_CLI
def test_cli_no_command_is_usage_error():
    assert _run().returncode == 2


@_CLI
def test_cli_missing_config_is_usage_error(tmp_path):
    r = _run("pf", "does_not_exist.xml", cwd=tmp_path)
    assert r.returncode == 2
    assert "no such config file" in r.stderr


@_CLI
@pytest.mark.integration
def test_cli_converged_exits_zero(tests_data_dir):
    r = _run("pf", "input_14.xml", "-q", "--no-timer", cwd=tests_data_dir)
    assert r.returncode == 0, r.stderr[-2000:]


@_CLI
@pytest.mark.integration
def test_cli_diverged_exits_three(tests_data_dir, tmp_path):
    """A diverged solve must not report success."""
    cfg = _capped(tests_data_dir, tmp_path)
    r = _run("pf", cfg.name, "-q", "--no-timer", cwd=tmp_path)
    assert r.returncode == 3, r.stderr[-2000:]
    assert "did not converge" in r.stderr


@_CLI
@pytest.mark.integration
def test_cli_honors_xml_nonlinear_false(tests_data_dir, tmp_path):
    """cursor.get() returns strings, so "false" must not select nl_solve.

    PETSc's SNES banner only appears on the nonlinear path, which is the
    only externally visible difference between the two solvers.
    """
    src = (tests_data_dir / "input_14.xml").read_text()
    cfg = tmp_path / "nlfalse.xml"
    cfg.write_text(src.replace("<tolerance>",
                               "<UseNonLinear>false</UseNonLinear>\n    <tolerance>"))
    shutil.copy(tests_data_dir / "IEEE14.raw", tmp_path / "IEEE14.raw")

    off = _run("pf", cfg.name, "--no-timer", cwd=tmp_path)
    assert off.returncode == 0, off.stderr[-2000:]
    assert "SNES" not in off.stdout + off.stderr

    on = _run("pf", "input_14.xml", "--solver", "nl", "--no-timer",
              cwd=tests_data_dir)
    assert "SNES" in on.stdout + on.stderr, "--solver nl no longer nonlinear"


@_CLI
@pytest.mark.integration
def test_cli_ca_without_contingencies_is_usage_error(tmp_path):
    """Analyzing nothing must not look like a clean run."""
    for f in (_CA_INPUT, _RAW_14):
        if not f.exists():
            pytest.skip(f"missing {f}")
    shutil.copy(_CA_INPUT, tmp_path / "input_14.xml")
    shutil.copy(_RAW_14, tmp_path / "IEEE14_ca.raw")
    r = _run("ca", "input_14.xml", "-q", "--no-timer", "--no-print-calcs",
             cwd=tmp_path)
    assert r.returncode == 2, r.stderr[-2000:]
    assert "FullBranchN1" in r.stderr


@_CLI
@pytest.mark.integration
@pytest.mark.mpi
def test_cli_diverged_exits_three_under_mpi(tests_data_dir, tmp_path,
                                            require_mpiexec):
    """Every rank must agree, and mpiexec must surface the code."""
    cfg = _capped(tests_data_dir, tmp_path)
    r = _run("pf", cfg.name, "-q", "--no-timer", cwd=tmp_path, np=2)
    assert r.returncode == 3, r.stderr[-2000:]
    assert r.stderr.count("did not converge") == 2


@_CLI
@pytest.mark.integration
@pytest.mark.mpi
def test_cli_dsf_finalizes_mpi(dsf_build_dir, require_mpiexec):
    """os._exit skipped MPI_Finalize, so mpiexec reported exit 1 on success."""
    r = _run("dsf", "input_9b3g.xml", "-q", "--no-timer",
             cwd=dsf_build_dir, np=2)
    assert r.returncode == 0, r.stderr[-2000:]
    assert "finalize" not in r.stderr.lower()


@_CLI
@pytest.mark.integration
def test_cli_ca_runs(ca_case):
    """A contingency that diverges is a finding, not a failed run."""
    r = _run("ca", "input_14.xml", "--no-timer", cwd=ca_case)
    assert r.returncode == 0, r.stderr[-2000:]
    assert "LINE_2_3" in r.stdout and "GEN_2" in r.stdout
    assert "DIVERGENT" in r.stdout


@_CLI
@pytest.mark.integration
@pytest.mark.mpi
def test_cli_ca_under_mpi_exits_zero(ca_case, require_mpiexec):
    """Regression: ranks solving different contingencies used to deadlock.

    Rank 1 reported "No reference bus found" and the run hung in solve().
    """
    r = _run("ca", "input_14.xml", "--no-timer", cwd=ca_case, np=2)
    assert r.returncode == 0, r.stderr[-2000:]
    assert "No reference bus found" not in r.stdout + r.stderr
    # gather() gives rank 0 every contingency, not just its own share.
    for name in ("LINE_2_3", "LINE_6_13", "LINE_13_14", "GEN_2"):
        assert name in r.stdout


@_CLI
@pytest.mark.integration
def test_cli_se_rejects_output_file(se_data_dir):
    """-o used to AttributeError after the solve had already printed."""
    r = _run("se", "input.xml", "-o", "se.out", cwd=se_data_dir)
    assert r.returncode == 2
    assert "no open()/close()" in r.stderr


@_CLI
def test_cli_rejects_unknown_psse_version(tests_data_dir, tmp_path):
    """An unknown format used to be ignored, exporting nothing and exiting 0."""
    r = _run("powerflow", str(tests_data_dir / "input_14.xml"),
             "--export-psse", "v99", str(tmp_path / "out.raw"))
    assert r.returncode == 2
    assert "--export-psse format must be one of" in r.stderr
    assert not (tmp_path / "out.raw").exists()
