# -------------------------------------------------------------
# file: python/tests/test_cli.py
# -------------------------------------------------------------
# CLI smoke tests: --help, --version.  A live PF run through the CLI
# is exercised by test_powerflow's integration case.
# -------------------------------------------------------------

import shutil
import subprocess

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
