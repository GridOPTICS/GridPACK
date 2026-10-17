# -------------------------------------------------------------
# file: python/tests/test_tools.py
# -------------------------------------------------------------
# gridpack.mcp.tools, directly and through the stdio server.
# -------------------------------------------------------------

import json
import shutil
import subprocess
import sys

import pytest

from gridpack.mcp.runner import JobRunner
from gridpack.mcp.tools import Tools

PF14 = "input/powerflow/input_14.xml"


@pytest.fixture(scope="module")
def tools(tmp_path_factory):
    return Tools(JobRunner(tmp_path_factory.mktemp("jobs")))


def test_list_cases(tools):
    cases = {c["case"]: c for c in tools.list_cases()}
    assert cases[PF14] == {"case": PF14, "apps": ["powerflow"],
                           "network": "IEEE14.raw"}
    assert cases["input/ca/input_14.xml"]["apps"] == ["powerflow", "ca"]
    # Unparsable, and a Dynamic_simulation-only input dsf aborts on.
    assert "input/ds/input_dsim.xml" not in cases
    assert "input/ds/input.xml" not in cases
    assert [c["case"] for c in tools.list_cases("se")] == [
        "input/se/input_118.xml", "input/se/input_14.xml"]
    with pytest.raises(ValueError):
        tools.list_cases("emt")


def test_describe_case(tools):
    d = tools.describe_case("input/ds/input_9b3g.xml")
    assert d["apps"] == ["powerflow", "dsf"]
    assert d["files"] == {"Powerflow.networkConfiguration": "9b3g.raw",
                          "Dynamic_simulation.generatorParameters": "9b3g.dyr"}
    assert d["settings"]["dsf"] == {"time_step": "0.01", "end_time": "1.0"}
    # The key solve() reads, not the NonlinearSolver one.
    assert tools.describe_case("input/ca/input_14.xml")["settings"][
        "powerflow"] == {"tolerance": "1.0e-6", "max_iterations": "50"}


def test_cases_with_missing_files_are_not_listed(tmp_path):
    (tmp_path / "input").mkdir()
    (tmp_path / "raw").mkdir()
    (tmp_path / "raw/here.raw").write_text("")
    for name, raw in (("ok", "here.raw"), ("broken", "gone.raw")):
        (tmp_path / ("input/%s.xml" % name)).write_text(
            "<Configuration><Powerflow><networkConfiguration>%s"
            "</networkConfiguration></Powerflow></Configuration>" % raw)
    assert [c["case"] for c in Tools(root=tmp_path).list_cases()] \
        == ["input/ok.xml"]


@pytest.mark.parametrize("kwargs", [
    {"case": "../../etc/passwd"},
    {"case": "/etc/passwd"},
    {"case": "input/se/input_14.xml"},           # se only
    {"case": PF14, "tolerance": -1.0},
    {"case": PF14, "ranks": 0},
])
def test_run_rejects_without_a_job(tools, kwargs):
    before = len(tools.runner.jobs())
    with pytest.raises(ValueError):
        tools.run_powerflow(**kwargs)
    assert len(tools.runner.jobs()) == before


@pytest.mark.integration
def test_run_powerflow_and_results(tools):
    from .test_powerflow import BUS5_VMAG_REF
    st = tools.run_powerflow(PF14)
    assert (st["state"], st["exit_code"], st["converged"]) == ("done", 0, True)
    assert "workdir" not in st
    res = tools.job_result(st["id"])
    bus5 = next(b for b in res["buses"] if b["busId"] == 5)
    assert abs(bus5["voltage"] - BUS5_VMAG_REF) < 1e-5 and "config" not in res
    assert tools.job_log(st["id"], "stdout", tail=3).count("\n") <= 2

    st = tools.run_powerflow(PF14, max_iterations=1)
    assert (st["state"], st["exit_code"], st["converged"]) == ("done", 3, False)


@pytest.mark.integration
def test_run_contingency_matches_ca_x(tools):
    from .conftest import CA_AUTO_N1_STATUS
    st = tools.run_contingency("input/ca/input_14.xml")
    res = tools.job_result(st["id"])
    assert {c["name"]: c["status"] for c in res["contingencies"]} \
        == CA_AUTO_N1_STATUS


@pytest.mark.integration
def test_run_dynamic_sim_and_state_estimation(tools):
    st = tools.run_dynamic_sim("input/ds/input_9b3g.xml", end_time=0.5)
    res = tools.job_result(st["id"])
    assert st["state"] == "done" and len(res["watch"]["rows"]) == 50

    st = tools.run_state_estimation("input/se/input_14.xml")
    res = tools.job_result(st["id"])
    assert st["converged"] is True and len(res["buses"]) == 14


def test_waiting_and_unknown_jobs(tmp_path):
    slow = tmp_path / "slow.py"
    slow.write_text("import time; time.sleep(30)")
    t = Tools(JobRunner(tmp_path / "jobs", cli=(sys.executable, str(slow))))
    st = t.run_powerflow(PF14, wait_seconds=0)
    assert st["state"] in ("queued", "running")
    with pytest.raises(ValueError, match="still"):
        t.job_result(st["id"])
    t.runner.cancel(st["id"])
    st = t.job_status(t.runner.wait(st["id"], 30).id)
    assert (st["state"], st["error"]) == ("failed", "cancelled")
    with pytest.raises(ValueError, match="no results"):
        t.job_result(st["id"])
    with pytest.raises(ValueError):
        t.job_log(st["id"], "env")
    with pytest.raises(KeyError):
        t.job_status("nope")


@pytest.mark.integration
@pytest.mark.skipif(shutil.which("gridpack-mcp") is None,
                    reason="gridpack-mcp console script not installed")
def test_tools_over_stdio():
    pytest.importorskip("mcp")
    p = subprocess.Popen(["gridpack-mcp"], stdin=subprocess.PIPE,
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                         text=True)

    def rpc(i, method, params):
        p.stdin.write(json.dumps({"jsonrpc": "2.0", "id": i,
                                  "method": method, "params": params}) + "\n")
        p.stdin.flush()
        return json.loads(p.stdout.readline())["result"]

    try:
        rpc(1, "initialize", {"protocolVersion": "2025-06-18",
                              "capabilities": {},
                              "clientInfo": {"name": "pytest", "version": "0"}})
        p.stdin.write(json.dumps({"jsonrpc": "2.0",
                                  "method": "notifications/initialized"}) + "\n")
        run = rpc(2, "tools/call", {"name": "run_powerflow",
                                    "arguments": {"case": PF14}})
        job = run["structuredContent"]["result"]
        res = rpc(3, "tools/call", {"name": "job_result",
                                    "arguments": {"job_id": job["id"]}})
        bad = rpc(4, "tools/call", {"name": "run_powerflow",
                                    "arguments": {"case": "../x"}})
        p.stdin.close()
        assert p.wait(timeout=30) == 0
    finally:
        p.kill()
    assert job["state"] == "done" and not run["isError"]
    assert len(res["structuredContent"]["result"]["buses"]) == 14
    # The reason reaches the client, not just "Error executing tool".
    assert bad["isError"] and "unknown case '../x'" in bad["content"][0]["text"]
