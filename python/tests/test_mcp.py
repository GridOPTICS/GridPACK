# -------------------------------------------------------------
# file: python/tests/test_mcp.py
# -------------------------------------------------------------
# gridpack-mcp: install, version and the stdio handshake.
# -------------------------------------------------------------

import json
import shutil
import subprocess
import sys

import pytest

_SCRIPT = shutil.which("gridpack-mcp")
_CLI = pytest.mark.skipif(_SCRIPT is None,
                          reason="gridpack-mcp console script not installed")


@_CLI
def test_version_is_the_package_version():
    from importlib.metadata import version
    r = subprocess.run([_SCRIPT, "--version"], capture_output=True,
                       text=True, timeout=60)
    assert r.returncode == 0
    assert r.stdout.strip() == "gridpack-mcp %s" % version("gridpack")


@_CLI
def test_stdio_handshake():
    pytest.importorskip("mcp")
    p = subprocess.Popen([_SCRIPT], stdin=subprocess.PIPE,
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                         text=True)

    def send(msg):
        p.stdin.write(json.dumps(msg) + "\n")
        p.stdin.flush()

    try:
        send({"jsonrpc": "2.0", "id": 1, "method": "initialize",
              "params": {"protocolVersion": "2025-06-18", "capabilities": {},
                         "clientInfo": {"name": "pytest", "version": "0"}}})
        # Any stray print on stdout breaks this parse.
        init = json.loads(p.stdout.readline())
        send({"jsonrpc": "2.0", "method": "notifications/initialized"})
        send({"jsonrpc": "2.0", "id": 2, "method": "tools/list"})
        tools = json.loads(p.stdout.readline())
        send({"jsonrpc": "2.0", "id": 3, "method": "prompts/list"})
        prompts = json.loads(p.stdout.readline())
        send({"jsonrpc": "2.0", "id": 4, "method": "prompts/get",
              "params": {"name": "n1_screen",
                         "arguments": {"case": "input/ds/input_9b3g.xml"}}})
        prompt = json.loads(p.stdout.readline())
        p.stdin.close()
        assert p.wait(timeout=30) == 0
        for line in p.stdout.read().splitlines():
            json.loads(line)
    finally:
        p.kill()
    assert init["result"]["serverInfo"]["name"] == "gridpack"
    assert [t["name"] for t in tools["result"]["tools"]] == [
        "list_cases", "describe_case", "create_case", "modify_case",
        "create_measurements", "run_powerflow", "run_contingency", "run_dynamic_sim",
        "run_state_estimation", "run_hadrec", "run_emt", "job_status",
        "job_result", "job_log", "job_series", "job_report", "compare_jobs",
        "plot_job", "run_study", "study_status", "study_result"]
    assert [q["name"] for q in prompts["result"]["prompts"]] == [
        "stability_study", "n1_screen", "voltage_study"]
    text = prompt["result"]["messages"][0]["content"]["text"]
    assert "run_study kind dynamic_n1 on input/ds/input_9b3g.xml" in text
    assert "max_ranks" in text and "replace_events" in text


@_CLI
@pytest.mark.parametrize("value", ["0", "-5", "nan", "soon"])
def test_timeout_must_be_positive(value):
    r = subprocess.run([_SCRIPT, "--timeout", value], capture_output=True,
                       text=True, timeout=60)
    assert r.returncode == 2 and "--timeout" in r.stderr


def test_no_cases_is_a_usage_error(monkeypatch, capsys):
    import gridpack.mcp.tools as mod
    from gridpack.mcp.server import main

    def missing():
        raise FileNotFoundError("no data_sets")
    monkeypatch.setattr(mod, "data_sets", missing)
    assert main([]) == 2
    assert "GRIDPACK_TEST_DATA" in capsys.readouterr().err


@_CLI
def test_missing_cases_directory_is_a_usage_error(tmp_path):
    r = subprocess.run([_SCRIPT, "--cases", str(tmp_path / "nope")],
                       capture_output=True, text=True, timeout=60)
    assert r.returncode == 2 and "does not exist" in r.stderr


@_CLI
@pytest.mark.parametrize("value", ["0", "-1", "two"])
def test_parallel_must_be_a_positive_integer(value):
    r = subprocess.run([_SCRIPT, "--parallel", value], capture_output=True,
                       text=True, timeout=60)
    assert r.returncode == 2 and "--parallel" in r.stderr


@_CLI
def test_doctor_checks_the_installation(tmp_path):
    r = subprocess.run([_SCRIPT, "doctor"], capture_output=True, text=True,
                       timeout=120)
    assert r.returncode == 0 and r.stdout.rstrip().endswith("ready")
    assert "ok       gridpack CLI" in r.stdout and "bundled cases" in r.stdout
    assert "mcp package" in r.stdout and "Python" in r.stdout
    r = subprocess.run([_SCRIPT, "doctor", "--cases", str(tmp_path / "nope")],
                       capture_output=True, text=True, timeout=120)
    assert r.returncode == 1 and "does not exist" in r.stdout
    assert r.stdout.rstrip().endswith("not ready")


@_CLI
def test_config_prints_both_registrations(tmp_path):
    r = subprocess.run([_SCRIPT, "config", "--parallel", "3", "--cases",
                        str(tmp_path)], capture_output=True, text=True,
                       timeout=60)
    assert r.returncode == 0
    head, js = r.stdout.split("# Claude Desktop", 1)
    assert head.startswith("# Claude Code\nclaude mcp add gridpack -e PATH=")
    assert "gridpack-mcp --parallel 3 --cases %s" % tmp_path in head
    cfg = json.loads(js.split("\n", 1)[1])["mcpServers"]["gridpack"]
    assert cfg["command"].endswith("gridpack-mcp")
    assert cfg["args"] == ["--parallel", "3", "--cases", str(tmp_path)]
    assert cfg["env"]["PATH"].split(":")[0] == str(
        __import__("pathlib").Path(_SCRIPT).parent)


def test_missing_mcp_is_a_usage_error():
    r = subprocess.run(
        [sys.executable, "-c",
         "import sys; sys.modules['mcp'] = None\n"
         "from gridpack.mcp.server import main; sys.exit(main([]))"],
        capture_output=True, text=True, timeout=60)
    assert r.returncode == 2
    assert "pip install 'gridpack[mcp]'" in r.stderr
