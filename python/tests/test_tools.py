# -------------------------------------------------------------
# file: python/tests/test_tools.py
# -------------------------------------------------------------
# gridpack.mcp.tools, directly and through the stdio server.
# -------------------------------------------------------------

import json
import os
import shutil
import subprocess
import sys
import time

import pytest

from gridpack.mcp.runner import JobRunner
from gridpack.mcp.tools import Tools

PF14 = "input/powerflow/input_14.xml"


@pytest.fixture(scope="module")
def tools(tmp_path_factory):
    t = Tools(JobRunner(tmp_path_factory.mktemp("jobs")))
    yield t
    # Jobs run in their own sessions and would outlive pytest.
    for job in t.runner.jobs():
        t.runner.cancel(job.id)


@pytest.fixture
def harmless(tmp_path):
    """Tools whose CLI and mpiexec are no-ops: a validation bug can't start
    a real run, such as a 9-bus case on 4 ranks."""
    noop = tmp_path / "noop.py"
    noop.write_text("")
    cmd = (sys.executable, str(noop))
    return Tools(JobRunner(tmp_path / "jobs", cli=cmd, mpiexec=cmd))


def test_list_cases(tools):
    cases = {c["case"]: c for c in tools.list_cases()}
    assert cases[PF14] == {"case": PF14, "apps": ["powerflow"],
                           "network": "IEEE14.raw", "buses": 14,
                           "max_ranks": min(os.cpu_count(), 14)}
    assert cases["input/ca/input_14.xml"]["apps"] == ["powerflow", "ca"]
    # Unparsable, and a Dynamic_simulation-only input dsf aborts on.
    assert "input/ds/input_dsim.xml" not in cases
    assert "input/ds/input.xml" not in cases
    assert [c["case"] for c in tools.list_cases("se")] == [
        "input/se/input_118.xml", "input/se/input_14.xml"]
    with pytest.raises(ValueError):
        tools.list_cases("kalman")


def test_hadrec_and_emt_cases(tools):
    cases = {c["case"]: c["apps"] for c in tools.list_cases()}
    # No fault event: HADREC would index an empty list.
    assert cases["input/ds/input_kundur_two_area_wshygp.xml"] == [
        "powerflow", "dsf"]
    assert cases["input/ds/input_9b3g.xml"] == ["powerflow", "dsf", "hadrec"]
    emt = [c["case"] for c in tools.list_cases("emt")]
    assert "emt/9-bus/input_9bus2.xml" in emt and len(emt) == 7
    assert tools.describe_case("emt/9-bus/input_9bus2.xml")["settings"][
        "emt"] == {"time_step": "0.00005", "end_time": "2.0"}


def test_describe_case(tools):
    d = tools.describe_case("input/ds/input_9b3g.xml")
    assert d["apps"] == ["powerflow", "dsf", "hadrec"]
    assert d["files"] == {"Powerflow.networkConfiguration": "9b3g.raw",
                          "Dynamic_simulation.generatorParameters": "9b3g.dyr"}
    assert d["settings"]["dsf"] == {"time_step": "0.01", "end_time": "1.0"}
    # The key solve() reads, not the NonlinearSolver one.
    assert tools.describe_case("input/ca/input_14.xml")["settings"][
        "powerflow"] == {"tolerance": "1.0e-6", "max_iterations": "50"}


def test_bus_counts_and_rank_caps(tools):
    cases = {c["case"]: c for c in tools.list_cases()}
    cpus = os.cpu_count()
    # Each count is the one the power flow reports for that network.
    for case, buses, cap in [
            ("emt/two-bus/input_2bus.xml", 2, 1),
            ("input/ds/input_9b3g.xml", 9, 3),
            ("input/ds/input_kundur_two_area_wsieg1.xml", 11, min(cpus, 11)),
            ("input/powerflow/input_240_v36.xml", 243, cpus),
            ("input/powerflow/input_european.xml", 13659, cpus)]:
        assert (cases[case]["buses"], cases[case]["max_ranks"]) \
            == (buses, cap), case


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
    {"case": "input/ds/input_kundur_two_area_wshygp.xml", "app": "hadrec"},
    # Over the cap.
    {"case": "input/ds/input_9b3g.xml", "ranks": 4},
    {"case": "emt/two-bus/input_2bus.xml", "ranks": 2},
])
def test_run_rejects_without_a_job(harmless, kwargs):
    run = getattr(harmless, {"hadrec": "run_hadrec"}.get(
        kwargs.pop("app", None), "run_powerflow"))
    with pytest.raises(ValueError):
        run(**kwargs)
    assert harmless.runner.jobs() == []


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


@pytest.mark.integration
def test_run_hadrec_and_emt(tools):
    st = tools.run_hadrec("input/ds/input_kundur_two_area_wsieg1.xml",
                          end_time=0.2)
    obs = tools.job_result(st["id"])["observations"]
    assert st["state"] == "done" and len(obs["columns"]) == 13
    assert abs(obs["rows"][-1][0] - 0.2) < 1e-9

    st = tools.run_emt("emt/9-bus/input_9bus2.xml", end_time=0.005)
    assert (st["state"], st["exit_code"]) == ("done", 0)
    job = tools.runner.get(st["id"])
    assert (job.workdir / ".petscrc").is_file()
    # EMT prints its timers even with --no-timer, so not the tail.
    assert "TS dt 5e-05 time 0.005" in job.log()


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

    t = Tools(JobRunner(tmp_path / "jobs2", cli=(sys.executable, str(slow))),
              timeout=1)
    st = t.run_powerflow(PF14)
    assert (st["state"], st["error"]) == ("failed", "timed out after 1 s")
    with pytest.raises(ValueError, match="no results"):
        t.job_result(st["id"])
    with pytest.raises(ValueError):
        t.job_log(st["id"], "env")
    with pytest.raises(KeyError):
        t.job_status("nope")


class _Server:
    """gridpack-mcp over stdio, initialised."""

    def __init__(self, *args):
        self.p = subprocess.Popen(["gridpack-mcp", *args], stdin=subprocess.PIPE,
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE, text=True)
        self.n = 0
        self.rpc("initialize", {"protocolVersion": "2025-06-18",
                                "capabilities": {},
                                "clientInfo": {"name": "pytest",
                                               "version": "0"}})
        self._send({"jsonrpc": "2.0", "method": "notifications/initialized"})

    def _send(self, msg):
        self.p.stdin.write(json.dumps(msg) + "\n")
        self.p.stdin.flush()

    def rpc(self, method, params):
        self.n += 1
        self._send({"jsonrpc": "2.0", "id": self.n, "method": method,
                    "params": params})
        return json.loads(self.p.stdout.readline())["result"]

    def call(self, name, **arguments):
        return self.rpc("tools/call", {"name": name, "arguments": arguments})

    def close(self):
        self.p.stdin.close()
        return self.p.wait(timeout=60)


_STDIO = pytest.mark.skipif(shutil.which("gridpack-mcp") is None,
                            reason="gridpack-mcp console script not installed")


def _procs(pattern):
    # pgrep -f matches the full command line; ps truncates it on Linux.
    out = subprocess.run(["pgrep", "-f", pattern], capture_output=True,
                         text=True).stdout
    return out.split()


@pytest.mark.integration
@_STDIO
def test_tools_over_stdio():
    pytest.importorskip("mcp")
    s = _Server()
    try:
        run = s.call("run_powerflow", case=PF14)
        job = run["structuredContent"]["result"]
        res = s.call("job_result", job_id=job["id"])
        bad = s.call("run_powerflow", case="../x")
        assert s.close() == 0
    finally:
        s.p.kill()
    assert job["state"] == "done" and not run["isError"]
    assert len(res["structuredContent"]["result"]["buses"]) == 14
    # The reason reaches the client, not just "Error executing tool".
    assert bad["isError"] and "unknown case '../x'" in bad["content"][0]["text"]


@pytest.mark.integration
@_STDIO
def test_server_exit_stops_its_jobs():
    pytest.importorskip("mcp")
    pattern = "dsf input_3000.xml --no-timer --json result.json"
    assert not _procs(pattern)
    s = _Server()
    try:
        st = s.call("run_dynamic_sim", case="input/ds/input_3000.xml",
                    end_time=1e5, wait_seconds=1)["structuredContent"]["result"]
        assert st["state"] == "running" and _procs(pattern)
        assert s.close() == 0
    finally:
        s.p.kill()
    for _ in range(50):
        if not _procs(pattern):
            break
        time.sleep(0.2)
    assert not _procs(pattern)


@pytest.mark.integration
@_STDIO
def test_server_timeout_kills_jobs():
    pytest.importorskip("mcp")
    s = _Server("--timeout", "2")
    try:
        st = s.call("run_dynamic_sim", case="input/ds/input_3000.xml",
                    end_time=1e5, wait_seconds=10)["structuredContent"]["result"]
        assert s.close() == 0
    finally:
        s.p.kill()
    assert (st["state"], st["error"]) == ("failed", "timed out after 2.0 s")


# -------------------------------------------------------------
# --cases directories
# -------------------------------------------------------------

def _models(tmp_path, name="mymodels"):
    from .conftest import _data_sets
    d, ds = tmp_path / name, _data_sets()
    (d / "sub").mkdir(parents=True)
    shutil.copy(ds / "input/powerflow/input_14.xml", d / "feeder.xml")
    shutil.copy(ds / "raw/IEEE14.raw", d)
    shutil.copy(ds / "input/ds/input_9b3g.xml", d / "sub/dyn.xml")
    for f in ("raw/9b3g.raw", "dyr/9b3g.dyr"):
        shutil.copy(ds / f, d / "sub")
    shutil.copy(ds / "contingencies/contingencies_14.xml", d)   # not a case
    (d / "broken.xml").write_text(
        "<Configuration><Powerflow><networkConfiguration>gone.raw"
        "</networkConfiguration></Powerflow></Configuration>")
    return d


def test_cases_directory_is_listed(tmp_path, capsys):
    t = Tools(JobRunner(tmp_path / "jobs"), cases=[_models(tmp_path)])
    mine = {c["case"]: c for c in t.list_cases()
            if c["case"].startswith("mymodels/")}
    assert sorted(mine) == ["mymodels/feeder.xml", "mymodels/sub/dyn.xml"]
    assert (mine["mymodels/sub/dyn.xml"]["buses"],
            mine["mymodels/sub/dyn.xml"]["max_ranks"]) == (9, 3)
    assert "_xml" not in mine["mymodels/feeder.xml"]
    err = capsys.readouterr().err
    assert "skipping mymodels/broken.xml: data files not found: gone.raw" \
        in err and "contingencies_14" not in err
    assert len([c for c in t.list_cases()
                if not c["case"].startswith("mymodels/")]) == 59


def test_cases_directory_rank_cap(tmp_path):
    noop = tmp_path / "noop.py"
    noop.write_text("")
    cmd = (sys.executable, str(noop))
    t = Tools(JobRunner(tmp_path / "jobs", cli=cmd, mpiexec=cmd),
              cases=[_models(tmp_path)])
    with pytest.raises(ValueError, match="at most 3 ranks"):
        t.run_dynamic_sim("mymodels/sub/dyn.xml", ranks=4)
    assert t.runner.jobs() == []


@pytest.mark.integration
def test_cases_directory_runs(tmp_path):
    from .test_powerflow import BUS5_VMAG_REF
    t = Tools(JobRunner(tmp_path / "jobs"), cases=[_models(tmp_path)])
    st = t.run_powerflow("mymodels/feeder.xml")
    bus5 = next(b for b in t.job_result(st["id"])["buses"] if b["busId"] == 5)
    assert st["converged"] and abs(bus5["voltage"] - BUS5_VMAG_REF) < 1e-5
    st = t.run_dynamic_sim("mymodels/sub/dyn.xml", end_time=0.2)
    assert (st["state"], st["exit_code"]) == ("done", 0)


def test_emt_petscrc_prefers_the_case_directory(tmp_path):
    from .conftest import _data_sets
    ds = _data_sets()
    for name, local in (("plain", False), ("own", True)):
        d = tmp_path / name
        shutil.copytree(ds / "emt/9-bus", d)
        if local:
            (d / ".petscrc").write_text("-emt_ts_type beuler\n")
        t = Tools(cases=[d])
        entry = t._entry("%s/input_9bus2.xml" % name)
        assert entry["_petscrc"] == (d.resolve() / ".petscrc" if local
                                     else ds / "emt/.petscrc")


@pytest.mark.parametrize("name, twice", [("input", False), ("emt", False),
                                         ("mymodels", True)])
def test_cases_directory_names_must_be_free(tmp_path, name, twice):
    d = _models(tmp_path, name)
    with pytest.raises(ValueError, match="is taken"):
        Tools(cases=[d, d] if twice else [d])
    with pytest.raises(ValueError, match="does not exist"):
        Tools(cases=[tmp_path / "nope"])


def test_without_data_sets_only_cases_directories(tmp_path, monkeypatch):
    import gridpack.mcp.tools as mod

    def missing():
        raise FileNotFoundError("no data_sets")
    monkeypatch.setattr(mod, "data_sets", missing)
    with pytest.raises(ValueError, match="GRIDPACK_TEST_DATA"):
        Tools()
    t = Tools(cases=[_models(tmp_path)])
    assert [c["case"] for c in t.list_cases()] == [
        "mymodels/feeder.xml", "mymodels/sub/dyn.xml"]
