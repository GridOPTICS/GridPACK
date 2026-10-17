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
    assert "<maxIteration>50</maxIteration>" in src, (
        "input_14.xml no longer has maxIteration=50; update this helper"
    )
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
    assert "<UseNonLinear>false</UseNonLinear>" in cfg.read_text(), (
        "config edit did nothing; the SNES assertion below would pass vacuously"
    )
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
    cfg = _CA_INPUT.read_text()
    for tag in ("FullBranchN1", "FullGeneratorN1"):
        on = "<%s>true</%s>" % (tag, tag)
        assert on in cfg, "CA input no longer sets %s" % tag
        cfg = cfg.replace(on, "<%s>false</%s>" % (tag, tag))
    (tmp_path / "input_14.xml").write_text(cfg)
    shutil.copy(_RAW_14, tmp_path / "IEEE14_ca.raw")
    r = _run("ca", "input_14.xml", "-q", "--no-timer", "--no-print-calcs",
             cwd=tmp_path)
    assert r.returncode == 2, r.stderr[-2000:]
    assert "no contingencies to run" in r.stderr


@_CLI
@pytest.mark.integration
def test_cli_ca_auto_n1_runs(ca_auto_case):
    """FullBranchN1/FullGeneratorN1 expand to ca.x's 27 cases."""
    r = _run("ca", "input_14.xml", "--no-timer", "--no-print-calcs",
             cwd=ca_auto_case)
    assert r.returncode == 0, r.stderr[-2000:]
    assert "Total contingencies: 27" in r.stdout
    assert "BR_13_14_B3" in r.stdout and "GN_8_1" in r.stdout


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
def test_cli_dsf_finalizes_mpi(dsf_data_dir, require_mpiexec):
    """os._exit skipped MPI_Finalize, so mpiexec reported exit 1 on success."""
    r = _run("dsf", "input_9b3g.xml", "-q", "--no-timer",
             cwd=dsf_data_dir, np=2)
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
    r = _run("se", "input_14.xml", "-o", "se.out", cwd=se_data_dir)
    assert r.returncode == 2
    assert "no open()/close()" in r.stderr


@_CLI
@pytest.mark.integration
def test_cli_dsf_quiet_suppresses_output(dsf_data_dir, tmp_path):
    """--quiet has to reach NoPrint; this XML sets no suppressOutput."""
    for name in ("9b3g.raw", "9b3g.dyr", "input_9b3g.xml"):
        shutil.copy(dsf_data_dir / name, tmp_path / name)
    loud = _run("dsf", "input_9b3g.xml", "--no-timer", cwd=tmp_path)
    quiet = _run("dsf", "input_9b3g.xml", "-q", "--no-timer", cwd=tmp_path)
    assert loud.returncode == 0 and quiet.returncode == 0, quiet.stderr[-2000:]
    assert "Monitoring generators" in loud.stdout
    assert "Monitoring generators" not in quiet.stdout


@_CLI
@pytest.mark.integration
def test_cli_hadrec_quiet_suppresses_output(dsf_data_dir, tmp_path):
    """--quiet has to reach NoPrint; teardown used to SEGV 139 here."""
    for name in ("9b3g.raw", "9b3g.dyr", "input_9b3g.xml"):
        shutil.copy(dsf_data_dir / name, tmp_path / name)
    loud = _run("hadrec", "input_9b3g.xml", "--no-timer", cwd=tmp_path)
    quiet = _run("hadrec", "input_9b3g.xml", "-q", "--no-timer", cwd=tmp_path)
    assert loud.returncode == 0 and quiet.returncode == 0, quiet.stderr[-2000:]
    assert "Power flow converged" in loud.stdout
    assert "Power flow converged" not in quiet.stdout
    assert "Monitoring generators" not in quiet.stdout


@_CLI
def test_cli_dsf_rejects_output_file(tmp_path):
    """-o wrote a 0-byte file and exited 0; open() captures nothing here."""
    cfg = tmp_path / "input.xml"
    cfg.write_text("<Configuration/>\n")
    r = _run("dsf", cfg.name, "-o", "dsf.out", cwd=tmp_path)
    assert r.returncode == 2
    assert "captures nothing" in r.stderr
    assert not (tmp_path / "dsf.out").exists()


@_CLI
def test_cli_rejects_unknown_psse_version(tests_data_dir, tmp_path):
    """An unknown format used to be ignored, exporting nothing and exiting 0."""
    r = _run("powerflow", str(tests_data_dir / "input_14.xml"),
             "--export-psse", "v99", str(tmp_path / "out.raw"))
    assert r.returncode == 2
    assert "--export-psse format must be one of" in r.stderr
    assert not (tmp_path / "out.raw").exists()


@_CLI
@pytest.mark.integration
def test_cli_emt_matches_emt_x(tmp_path):
    """Regression: on Linux the EMT control blocks bound to the dynamic-
    simulation module's same-named classes and the run diverged from step 1.
    """
    from gridpack._libpaths import EMT_MODULE
    emt_x = Path(EMT_MODULE).resolve().parent.parent / "bin" / "emt.x"
    if not emt_x.exists():
        pytest.skip(f"missing {emt_x}")
    src = _DATA / "emt"
    for d in ("driver", "cli"):
        run = tmp_path / d
        run.mkdir()
        for f in ("9-bus/input_9bus2.xml", "9-bus/case9_PSCAD.raw",
                  "9-bus/case9_GAST.dyr", ".petscrc"):
            shutil.copy(src / f, run / Path(f).name)
        cfg = run / "input_9bus2.xml"
        text = cfg.read_text()
        assert "<simulationTime>2.0</simulationTime>" in text
        # 200 steps: the corruption showed from the second row.
        cfg.write_text(text.replace("<simulationTime>2.0</simulationTime>",
                                    "<simulationTime>0.01</simulationTime>"))

    d = subprocess.run([str(emt_x), "input_9bus2.xml"], cwd=tmp_path / "driver",
                       capture_output=True, text=True, timeout=300)
    assert d.returncode == 0, d.stderr[-2000:]
    r = _run("emt", "input_9bus2.xml", cwd=tmp_path / "cli")
    assert r.returncode == 0, r.stderr[-2000:]
    assert ((tmp_path / "cli/emtoutput.csv").read_bytes()
            == (tmp_path / "driver/emtoutput.csv").read_bytes())
