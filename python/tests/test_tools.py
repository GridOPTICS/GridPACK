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
import xml.etree.ElementTree as ET

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
    assert all(c["notes"] == ["emt results are status only"]
               for c in tools.list_cases("emt"))
    assert "notes" not in tools.describe_case(PF14) and "notes" not in next(
        c for c in tools.list_cases("powerflow") if "emt" not in c["apps"])
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


# -------------------------------------------------------------
# create_case
# -------------------------------------------------------------

def _mine(tmp_path, *files):
    from .conftest import _data_sets
    d = tmp_path / "mine"
    (d / "dyn").mkdir(parents=True)
    for f in files:
        shutil.copy(_data_sets() / f, d / ("dyn" if "9b3g" in f else ""))
    return d


def test_create_case_from_raw_and_dyr(tmp_path, capsys):
    d = _mine(tmp_path, "raw/IEEE14.raw", "raw/9b3g.raw", "dyr/9b3g.dyr")
    t = Tools(JobRunner(tmp_path / "jobs"), cases=[d])
    pf = t.create_case("mine/IEEE14.raw")
    assert (pf["case"], pf["apps"], pf["buses"]) == (
        "mine/IEEE14.xml", ["powerflow"], 14)
    assert pf["files"] == {"Powerflow.networkConfiguration": "IEEE14.raw"}
    dyn = t.create_case(str(d / "dyn/9b3g.raw"), dyr="mine/dyn/9b3g.dyr",
                        name="nine")
    # hadrec needs a fault event, which the template has not got.
    assert (dyn["case"], dyn["apps"]) == ("mine/dyn/nine.xml",
                                          ["powerflow", "dsf"])
    assert (dyn["buses"], dyn["max_ranks"]) == (9, 3)
    assert dyn["files"]["Dynamic_simulation.generatorParameters"] == "9b3g.dyr"
    assert dyn["settings"]["dsf"] == {"time_step": "0.01", "end_time": "10.0"}
    assert {c["case"] for c in t.list_cases()} >= {"mine/IEEE14.xml",
                                                   "mine/dyn/nine.xml"}
    watch = ET.parse(d / "dyn/nine.xml").getroot().find(
        "Dynamic_simulation/generatorWatch")
    assert [(g.find("busID").text, g.find("generatorID").text)
            for g in watch] == [("1", "1"), ("2", "1"), ("3", "1")]
    with pytest.raises(ValueError, match="mine/IEEE14.xml exists"):
        t.create_case("mine/IEEE14.raw")
    assert capsys.readouterr().err == ""


@pytest.mark.parametrize("raw, tag, gens", [
    ("9b3g.raw", "networkConfiguration", 3),
    ("kundur-twoarea_v33.raw", "networkConfiguration_v33", 4),
    ("240busWECC_2018_PSS_fixedshunt.raw", "networkConfiguration_v34", 146),
    ("240busWECC_2018_PSS_PQLoad_fixedshunt_noremotebus_yuan_v36.raw",
     "networkConfiguration_v36", 146),
])
def test_raw_version_tag_and_generators(raw, tag, gens):
    from .conftest import _data_sets
    from gridpack.mcp.templates import build_xml, raw_generators, raw_version
    path = _data_sets() / "raw" / raw
    g = raw_generators(path)
    assert len(g) == gens and all(isinstance(b, int) and i for b, i in g)
    root = ET.fromstring(build_xml(raw, raw_version(path), dyr="x.dyr",
                                   generators=g))
    assert root.find("Powerflow/%s" % tag).text == raw
    assert len(root.find("Dynamic_simulation/generatorWatch")) == gens
    assert root.find("Dynamic_simulation/Events") is not None


def test_unlabelled_v33_generator_block(tmp_path):
    """PSS/E 33 labels its block ends; hand-made files may not."""
    from gridpack.mcp.templates import raw_generators, raw_version
    raw = tmp_path / "plain33.raw"
    raw.write_text("0, 100.00, 33, 0, 1, 60.00\ntitle\ntitle\n"
                   " 1,'A',1\n 2,'B',1\n0\n"          # buses
                   " 1,'1 ',1\n0\n"                    # loads
                   " 2,'1 ',1\n0\n"                    # fixed shunts
                   " 2,'G1',1\n 1,'2 ',1\n0\n"        # generators
                   " 1,2,'1 '\n0\nQ\n")
    assert raw_version(raw) == 33
    assert raw_generators(raw) == [(2, "G1"), (1, "2")]


def test_create_case_rejections(tmp_path):
    d = _mine(tmp_path, "raw/IEEE14.raw")
    outside = tmp_path / "IEEE14.raw"
    shutil.copy(d / "IEEE14.raw", outside)
    (d / "nogen.raw").write_text("0 100.0\n\n\n 1,'B',1\n0\n0\n")
    (d / "nogen.dyr").write_text("")
    t = Tools(JobRunner(tmp_path / "jobs"), cases=[d])
    for kwargs, msg in [
            ({"raw": str(outside)}, "RAW file must be under a cases dir"),
            ({"raw": "mine/../IEEE14.raw"}, "RAW file must be under"),
            ({"raw": "mine/gone.raw"}, "RAW file mine/gone.raw not found"),
            ({"raw": "mine/IEEE14.raw", "dyr": "mine/gone.dyr"},
             "DYR file mine/gone.dyr not found"),
            ({"raw": "mine/IEEE14.raw", "name": "../x"}, "name must be"),
            ({"raw": "mine/IEEE14.raw", "name": "a b"}, "name must be"),
            ({"raw": "mine/nogen.raw", "dyr": "mine/nogen.dyr"},
             "no generator records")]:
        with pytest.raises(ValueError, match=msg):
            t.create_case(**kwargs)
    assert not list(tmp_path.rglob("*.xml"))                # nothing written
    with pytest.raises(ValueError, match="--cases DIR"):
        Tools(JobRunner(tmp_path / "jobs2")).create_case("x.raw")


@pytest.mark.integration
def test_created_cases_run(tmp_path):
    from .test_powerflow import BUS5_VMAG_REF
    d = _mine(tmp_path, "raw/IEEE14.raw", "raw/9b3g.raw", "dyr/9b3g.dyr")
    t = Tools(JobRunner(tmp_path / "jobs"), cases=[d])
    st = t.run_powerflow(t.create_case("mine/IEEE14.raw")["case"])
    bus5 = next(b for b in t.job_result(st["id"])["buses"] if b["busId"] == 5)
    assert st["converged"] and abs(bus5["voltage"] - BUS5_VMAG_REF) < 1e-5
    dyn = t.create_case("mine/dyn/9b3g.raw", dyr="mine/dyn/9b3g.dyr")["case"]
    st = t.run_dynamic_sim(dyn, end_time=0.1)
    w = t.job_result(st["id"])["watch"]
    assert (st["state"], st["exit_code"], len(w["rows"])) == ("done", 0, 10)
    assert {c for c in w["columns"] if c.endswith("_speed")} == {
        "1_1_speed", "2_1_speed", "3_1_speed"}
    # Without events the system holds its power-flow point.
    i = w["columns"].index("1_1_speed")
    assert all(abs(r[i] - 1.0) < 1e-6 for r in w["rows"])


# -------------------------------------------------------------
# modify_case
# -------------------------------------------------------------

def test_modify_case_derives_and_lists(tmp_path):
    t = Tools(JobRunner(tmp_path / "jobs"))
    edits = [{"op": "scale_load", "bus": 9, "factor": 1.2},
             {"op": "set_branch_in_service", "from": 1, "to": 5, "ckt": "BL",
              "in_service": False}]
    d = t.modify_case(PF14, edits)
    assert d["case"] == "derived/input_14-1/input_14-1.xml"
    assert (d["base"], d["apps"], d["buses"]) == (PF14, ["powerflow"], 14)
    assert d["files"] == {"Powerflow.networkConfiguration": "IEEE14.raw"}
    assert [(e["was"], e["now"]) for e in d["edits"]] == [
        ({"PL": 29.5, "QL": 16.6}, {"PL": 35.4, "QL": 19.92}),
        ({"STAT": 1}, {"STAT": 0})]
    root = tmp_path / "jobs/derived/input_14-1"
    assert sorted(f.name for f in root.iterdir()) == [
        "IEEE14.raw", "edits.json", "input_14-1.xml"]
    assert json.loads((root / "edits.json").read_text())["base"] == PF14
    # The base case's file is untouched; the copy differs on two lines.
    from .conftest import _data_sets
    a = (_data_sets() / "raw/IEEE14.raw").read_text().splitlines()
    b = (root / "IEEE14.raw").read_text().splitlines()
    assert sum(x != y for x, y in zip(a, b)) == 2
    assert "derived/input_14-1/input_14-1.xml" in [
        c["case"] for c in t.list_cases()]
    # Names count up; a derived case can be derived again.
    d2 = t.modify_case(d["case"], [{"op": "set_generator_active_power",
                                    "bus": 2, "mw": 60.0}], name="again")
    assert d2["case"] == "derived/again/again.xml"
    assert d2["edits"][0]["was"] == {"PG": 40.0}
    assert t.modify_case(PF14, edits)["case"].startswith("derived/input_14-2/")


def test_modify_case_keeps_dynamics_and_caps(tmp_path):
    t = Tools(JobRunner(tmp_path / "jobs"))
    d = t.modify_case("input/ds/input_9b3g.xml", [
        {"op": "set_load_active_power", "bus": 5, "mw": 100.0}])
    assert (d["apps"], d["max_ranks"]) == (["powerflow", "dsf", "hadrec"], 3)
    assert d["files"] == {"Powerflow.networkConfiguration": "9b3g.raw",
                          "Dynamic_simulation.generatorParameters": "9b3g.dyr"}
    assert d["settings"]["dsf"] == {"time_step": "0.01", "end_time": "1.0"}


@pytest.mark.parametrize("kwargs, msg", [
    ({"case": "nope.xml", "edits": [{"op": "scale_load", "bus": 9,
                                     "factor": 2}]}, "unknown case"),
    ({"case": PF14, "edits": []}, "non-empty list"),
    ({"case": PF14, "edits": [{"op": "scale_load", "bus": 99, "factor": 2}]},
     "no bus 99"),
    ({"case": PF14, "edits": [{"op": "scale_load", "bus": 9, "factor": 2}],
     "name": "a/b"}, "name must be"),
])
def test_modify_case_rejections_leave_nothing(tmp_path, kwargs, msg):
    t = Tools(JobRunner(tmp_path / "jobs"))
    with pytest.raises(ValueError, match=msg):
        t.modify_case(**kwargs)
    assert not (tmp_path / "jobs/derived").exists() or \
        not list((tmp_path / "jobs/derived").iterdir())
    assert not any(c["case"].startswith("derived/") for c in t.list_cases())


def test_modify_case_edits_a_state_estimation_network(tmp_path):
    t = Tools(JobRunner(tmp_path / "jobs"))
    d = t.modify_case("input/se/input_14.xml", [
        {"op": "scale_load", "bus": 9, "factor": 2}])
    assert d["apps"] == ["se"] and d["edits"][0]["now"] == {"PL": 59.0,
                                                            "QL": 33.2}
    assert "State_estimation.networkConfiguration" in d["files"]


def test_derived_is_a_reserved_cases_name(tmp_path):
    with pytest.raises(ValueError, match="is taken"):
        Tools(cases=[_models(tmp_path, "derived")])


@pytest.mark.integration
def test_modified_cases_run(tmp_path):
    from .test_powerflow import BUS5_VMAG_REF
    t = Tools(JobRunner(tmp_path / "jobs"))
    d = t.modify_case(PF14, [
        {"op": "set_load_active_power", "bus": 9, "mw": 60.0},
        {"op": "set_generator_active_power", "bus": 2, "mw": 60.0},
        {"op": "set_branch_in_service", "from": 1, "to": 5, "ckt": "BL",
         "in_service": False},
        {"op": "set_shunt", "bus": 9, "mvar": 25.0}])
    st = t.run_powerflow(d["case"])
    res = t.job_result(st["id"])
    bus = {b["busId"]: b for b in res["buses"]}
    assert st["converged"] and bus[9]["pLoad"] == 60.0 and bus[9]["shuntMvar"] == 25.0
    assert next(g for g in res["generators"] if g["busId"] == 2)["pGen"] == 60.0
    assert not any({b["fromBus"], b["toBus"]} == {1, 5} for b in res["branches"])
    assert abs(bus[5]["voltage"] - BUS5_VMAG_REF) > 1e-3     # really changed

    # v33, and the dynamics still run on the derived case.
    k = t.modify_case("input/ds/input_kundur_two_area_esst1a.xml", [
        {"op": "set_load_active_power", "bus": 7, "mw": 900.0},
        {"op": "set_shunt_in_service", "bus": 9, "in_service": False}])
    st = t.run_powerflow(k["case"])
    bus = {b["busId"]: b for b in t.job_result(st["id"])["buses"]}
    assert st["converged"] and bus[7]["pLoad"] == 900.0 and bus[9]["shuntMvar"] == 0.0
    st = t.run_dynamic_sim(k["case"], end_time=0.05)
    assert (st["state"], st["exit_code"]) == ("done", 0)

    # v34 with labelled blocks and a lettered machine id.
    w = t.modify_case("input/input_240bus.xml", [
        {"op": "set_generator_in_service", "bus": 1032, "id": "C",
         "in_service": False}])
    st = t.run_powerflow(w["case"])
    gen = next(g for g in t.job_result(st["id"])["generators"]
               if g["busId"] == 1032 and g["genId"].strip() == "C")
    assert st["converged"] and gen["status"] == 0


# -------------------------------------------------------------
# modify_case: events
# -------------------------------------------------------------

DS9 = "input/ds/input_9b3g.xml"


def _events(xml):
    root = ET.parse(xml).getroot()
    return [(e.tag, {c.tag: c.text.strip() for c in e})
            for e in root.find("Dynamic_simulation/Events")]


def test_modify_case_adds_events_in_order(tmp_path):
    t = Tools(JobRunner(tmp_path / "jobs"))
    d = t.modify_case(DS9, [
        {"op": "add_fault", "from": 6, "to": 7, "start": 0.1, "duration": 0.05},
        {"op": "scale_load", "bus": 5, "factor": 1.1},
        {"op": "trip_branch", "from": 8, "to": 7, "time": 0.2},
        {"op": "trip_generator", "bus": 3, "time": 0.3},
        {"op": "add_bus_fault", "bus": 5, "start": 0.4, "duration": 0.05,
         "b_pu": 1000}])
    assert d["apps"] == ["powerflow", "dsf", "hadrec"]
    # The base case's own fault stays first; ids and the line's
    # orientation come from the RAW.
    assert _events(tmp_path / "jobs/derived/input_9b3g-1/input_9b3g-1.xml") == [
        ("faultEvent", {"beginFault": "0.03", "endFault": "0.06",
                        "faultBranch": "6 7", "timeStep": "0.01"}),
        ("faultEvent", {"beginFault": "0.1", "endFault": "0.15",
                        "faultBranch": "6 7", "id": "BL"}),
        ("LineStatus", {"time": "0.2", "line": "7 8", "id": "BL",
                        "status": "0"}),
        ("GenStatus", {"time": "0.3", "bus": "3", "id": "1", "status": "0"}),
        ("BusFault", {"begin": "0.4", "end": "0.45", "bus": "5",
                      "yfault": "0 1000"})]
    # The report keeps the input order, events as XML.
    assert [e["op"] for e in d["edits"]] == [
        "add_fault", "scale_load", "trip_branch", "trip_generator",
        "add_bus_fault"]
    assert d["edits"][1]["now"] == {"PL": 99.0, "QL": 33.0}
    assert [e["xml"].split(">")[0] for e in d["edits"] if "xml" in e] == [
        "<faultEvent", "<LineStatus", "<GenStatus", "<BusFault"]
    # A derived case with no network edit keeps the RAW byte for byte.
    e = t.modify_case(DS9, [{"op": "trip_generator", "bus": 2, "time": 0.5,
                             "in_service": True}], name="ev")
    from .conftest import _data_sets
    assert (tmp_path / "jobs/derived/ev/9b3g.raw").read_bytes() == \
        (_data_sets() / "raw/9b3g.raw").read_bytes()
    assert _events(tmp_path / "jobs/derived/ev/ev.xml")[-1] == (
        "GenStatus", {"time": "0.5", "bus": "2", "id": "1", "status": "1"})


def test_hadrec_needs_a_fault_first(tmp_path):
    d = _mine(tmp_path, "raw/9b3g.raw", "dyr/9b3g.dyr")
    t = Tools(JobRunner(tmp_path / "jobs"), cases=[d])
    plain = t.create_case("mine/dyn/9b3g.raw", dyr="mine/dyn/9b3g.dyr")["case"]
    trip = t.modify_case(plain, [{"op": "trip_generator", "bus": 3, "time": 1.0}])
    assert trip["apps"] == ["powerflow", "dsf"]
    fault = t.modify_case(plain, [
        {"op": "add_bus_fault", "bus": 7, "start": 1.0, "duration": 0.1},
        {"op": "trip_generator", "bus": 3, "time": 1.1}])
    assert fault["apps"] == ["powerflow", "dsf", "hadrec"]
    with pytest.raises(ValueError, match="cannot run hadrec"):
        t.run_hadrec(trip["case"])          # refused before any job starts
    assert t.runner.jobs() == []


def test_transformers_can_be_tripped_as_events(tmp_path):
    t = Tools(JobRunner(tmp_path / "jobs"))
    k = t.modify_case("input/ds/input_kundur_two_area_esst1a.xml", [
        {"op": "trip_branch", "from": 1, "to": 5, "time": 0.5}])
    assert _events(tmp_path / "jobs/derived" / k["case"].split("/", 1)[1])[-1] == (
        "LineStatus", {"time": "0.5", "line": "5 1", "id": "1", "status": "0"})


@pytest.mark.parametrize("case, edit, msg", [
    (DS9, {"op": "trip_generator", "bus": 3, "time": 0.105}, "not a multiple"),
    (DS9, {"op": "add_fault", "from": 6, "to": 7, "start": 0.98,
           "duration": 0.05}, "past simulationTime 1"),
    (DS9, {"op": "add_fault", "from": 6, "to": 7, "start": 0.1,
           "duration": 0}, "duration must be positive"),
    (DS9, {"op": "trip_generator", "bus": 3, "time": -0.1}, "not be negative"),
    (DS9, {"op": "add_bus_fault", "bus": 50, "start": 0.1, "duration": 0.1},
     "no bus 50"),
    (DS9, {"op": "trip_branch", "from": 1, "to": 9, "time": 0.1},
     "no branch 1-9"),
    (DS9, {"op": "add_fault", "from": 6, "to": 7, "ckt": "1", "start": 0.1,
           "duration": 0.1}, "has circuit BL, not 1"),
    (DS9, {"op": "trip_generator", "bus": 5, "time": 0.1},
     "no generator at bus 5"),
    (DS9, {"op": "trip_generator", "bus": 3, "id": "2", "time": 0.1},
     "generator id 1, not 2"),
    (DS9, {"op": "trip_branch", "from": 6, "to": 7, "time": 0.1,
           "in_service": 0}, "true or false"),
    (DS9, {"op": "trip_branch", "from": 6, "to": 7, "time": "0.1"},
     "time must be a number"),
    (DS9, {"op": "add_fault", "from": 6, "to": 7, "start": 0.1}, "missing duration"),
    (DS9, {"op": "shed_load", "bus": 5}, "unknown op 'shed_load'; one of .*add_fault"),
    (PF14, {"op": "add_fault", "from": 1, "to": 5, "start": 0.1,
            "duration": 0.1}, "no Dynamic_simulation block"),
    ("input/ds/input_kundur_two_area_esst1a.xml",
     {"op": "trip_branch", "from": 7, "to": 8, "time": 0.5},
     "circuits 1, 2; give ckt"),
])
def test_rejected_events_leave_nothing(tmp_path, case, edit, msg):
    t = Tools(JobRunner(tmp_path / "jobs"))
    with pytest.raises(ValueError, match=msg):
        t.modify_case(case, [{"op": "scale_load", "bus": 5, "factor": 1.1}
                             if case == DS9 else
                             {"op": "scale_load", "bus": 9, "factor": 1.1},
                             edit])
    assert not (tmp_path / "jobs/derived").exists()


@pytest.mark.integration
def test_events_change_the_run_and_hadrec_accepts_a_fault(tmp_path):
    t = Tools(JobRunner(tmp_path / "jobs"))
    base = t.job_result(t.run_dynamic_sim(DS9, end_time=0.3)["id"])["watch"]
    d = t.modify_case(DS9, [{"op": "trip_generator", "bus": 3, "time": 0.2}])
    st = t.run_dynamic_sim(d["case"], end_time=0.3)
    w = t.job_result(st["id"])["watch"]
    i = w["columns"].index("1_1_speed")
    same = [abs(a[i] - b[i]) for a, b in zip(w["rows"], base["rows"])]
    assert st["exit_code"] == 0 and len(same) == 30
    assert max(same[:19]) == 0 and max(same[21:]) > 1e-5    # trip at 0.2 s
    # Line and bus-fault events change the trajectory too.
    for edit in ({"op": "trip_branch", "from": 7, "to": 8, "time": 0.1},
                 {"op": "add_bus_fault", "bus": 5, "start": 0.1,
                  "duration": 0.05}):
        e = t.modify_case(DS9, [edit])
        w = t.job_result(t.run_dynamic_sim(e["case"], end_time=0.3)["id"])["watch"]
        assert max(abs(a[i] - b[i]) for a, b in zip(w["rows"], base["rows"])) > 1e-5
    # hadrec runs a derived case whose first event is a fault.
    mine = _mine(tmp_path, "raw/9b3g.raw", "dyr/9b3g.dyr")
    t2 = Tools(JobRunner(tmp_path / "jobs2"), cases=[mine])
    plain = t2.create_case("mine/dyn/9b3g.raw", dyr="mine/dyn/9b3g.dyr")["case"]
    f = t2.modify_case(plain, [{"op": "add_fault", "from": 6, "to": 7,
                                "start": 0.1, "duration": 0.05}])
    st = t2.run_hadrec(f["case"], end_time=0.3)
    obs = t2.job_result(st["id"])["observations"]
    assert (st["state"], st["exit_code"]) == ("done", 0) and len(obs["rows"]) > 0


# -------------------------------------------------------------
# post-processing tools
# -------------------------------------------------------------

def test_post_tools_need_a_finished_job(tmp_path):
    slow = tmp_path / "slow.py"
    slow.write_text("import time; time.sleep(30)")
    t = Tools(JobRunner(tmp_path / "jobs", cli=(sys.executable, str(slow))))
    st = t.run_powerflow(PF14, wait_seconds=0)
    for call in (lambda: t.job_series(st["id"]), lambda: t.job_report(st["id"]),
                 lambda: t.compare_jobs(st["id"], st["id"]),
                 lambda: t.plot_job(st["id"])):
        with pytest.raises(ValueError, match="still"):
            call()
    t.runner.cancel(st["id"])
    with pytest.raises(KeyError):
        t.job_report("nope")


@pytest.mark.integration
def test_post_tools_on_real_jobs(tmp_path):
    t = Tools(JobRunner(tmp_path / "jobs"))
    pf = t.run_powerflow(PF14)["id"]
    text = t.job_report(pf)
    assert text.startswith("Power flow converged in 2 iterations; 14 buses, 20 branches, 5 generators.")
    assert "Slack bus 1 supplies 232.4 MW" in text
    assert "No branch ratings in the case" in text
    opened = t.modify_case(PF14, [{"op": "set_branch_in_service", "from": 1,
                                   "to": 5, "ckt": "BL", "in_service": False}])
    pf2 = t.run_powerflow(opened["case"])["id"]
    c = t.compare_jobs(pf, pf2)
    assert c["branches"] == {"only_in_a": ["1-5 BL"], "only_in_b": []}
    assert c["voltage_pu"]["at_bus"] == 5 and c["voltage_pu"]["max_change"] < -0.01

    base = t.run_dynamic_sim(DS9, end_time=0.5)["id"]
    s = t.job_series(base, ["*_speed"], max_points=6)
    assert s["columns"] == ["t", "1_1_speed"] and s["points"] == 6
    assert s["decimated_from"] == 50 and s["rows"][-1][0] == pytest.approx(0.49)
    text = t.job_report(base)
    assert "Event: fault at bus 6 from 0.03 to 0.06 s (30 ms)" in text
    assert "Peak speed deviation" in text and "Verdict:" in text
    trip = t.modify_case(DS9, [{"op": "trip_generator", "bus": 3, "time": 0.2}])
    d = t.compare_jobs(base, t.run_dynamic_sim(trip["case"], end_time=0.5)["id"])
    assert d["first_difference_at"] == pytest.approx(0.2) and d["rows_compared"] == 50
    with pytest.raises(ValueError, match="cannot compare"):
        t.compare_jobs(pf, base)
    try:
        import matplotlib  # noqa: F401
    except ImportError:
        with pytest.raises(ValueError, match="gridpack\\[plot\\]"):
            t.plot_job(base)
    else:
        assert t.plot_job(base)[:4] == b"\x89PNG"


# -------------------------------------------------------------
# studies
# -------------------------------------------------------------

def test_modify_case_can_replace_events(tmp_path):
    t = Tools(JobRunner(tmp_path / "jobs"))
    d = t.modify_case(DS9, [{"op": "add_bus_fault", "bus": 7, "start": 0.2,
                             "duration": 0.05}], replace_events=True)
    xml = tmp_path / "jobs/derived" / d["case"].split("/", 1)[1]
    assert [tag for tag, _ in _events(xml)] == ["BusFault"]
    e = t.modify_case(DS9, [{"op": "scale_load", "bus": 5, "factor": 1.1}],
                      replace_events=True)
    assert _events(tmp_path / "jobs/derived" / e["case"].split("/", 1)[1]) == []
    assert e["apps"] == ["powerflow", "dsf"]           # no fault left for hadrec


def test_run_study_fans_out_before_running(tmp_path):
    noop = tmp_path / "noop.py"
    noop.write_text("")
    t = Tools(JobRunner(tmp_path / "jobs", max_parallel=2,
                        cli=(sys.executable, str(noop)),
                        mpiexec=(sys.executable, str(noop))))
    st = t.run_study(DS9, "dynamic_n1", end_time=0.5, wait_seconds=5)
    assert st["kind"] == "dynamic_n1" and st["items"] == 9      # 9b3g branches
    assert st["params"]["branches"][:2] == ["1-4 BL", "2-8 BL"]
    assert [j["state"] for j in st["jobs"]] == ["failed"] * 9    # no-op CLI
    derived = sorted(p.name for p in (tmp_path / "jobs/derived").iterdir())
    assert len(derived) == 9 and all(p.startswith(st["id"] + "-out-") for p in derived)
    xml = tmp_path / "jobs/derived" / st["jobs"][0]["case"].split("/", 1)[1]
    assert _events(xml) == [
        ("faultEvent", {"beginFault": "0.1", "endFault": "0.15",
                        "faultBranch": "1 4", "id": "BL"}),
        ("LineStatus", {"time": "0.15", "line": "1 4", "id": "BL", "status": "0"})]
    # Every machine is watched, whatever the base case watched.
    watch = ET.parse(xml).getroot().find("Dynamic_simulation/generatorWatch")
    assert [g.find("busID").text for g in watch] == ["1", "2", "3"]
    res = t.study_result(st["id"])
    assert (res["finished"], res["pending"], len(res["failed"])) == (0, 0, 9)
    assert t.study_status(st["id"])["states"] == {"failed": 9}
    sw = t.run_study(PF14, "load_sweep", factors=[1.0, 1.1], scale_generation=True,
                     wait_seconds=5)
    assert sw["items"] == 2 and sw["jobs"][1]["name"] == "load-1.1"
    edits = json.loads((tmp_path / "jobs/derived" / sw["jobs"][1]["case"].split("/", 1)[1]
                        ).with_name("edits.json").read_text())["edits"]
    assert [e["op"] for e in edits] == ["scale_all_loads", "scale_all_generation"]
    assert t.run_study(DS9, "dynamic_n1", branches=["4-5", "6-7 BL"], end_time=0.5,
                       wait_seconds=0)["items"] == 2
    with pytest.raises(KeyError):
        t.study_status("nope")


@pytest.mark.parametrize("kwargs, msg", [
    ({"case": DS9, "kind": "n-2"}, "kind is one of"),
    ({"case": PF14, "kind": "dynamic_n1"}, "cannot run dsf"),
    ({"case": DS9, "kind": "load_sweep", "factors": []}, "non-empty list"),
    ({"case": DS9, "kind": "load_sweep", "factors": [1.0, -2]}, "positive numbers"),
    ({"case": DS9, "kind": "load_sweep", "factors": "1.1"}, "non-empty list"),
    ({"case": DS9, "kind": "dynamic_n1", "ranks": 4}, "at most 3 ranks"),
    ({"case": DS9, "kind": "dynamic_n1", "ranks": 0}, "positive int"),
    ({"case": DS9, "kind": "dynamic_n1", "fault_start": 0.9, "end_time": 0.9},
     "must end before end_time"),
    ({"case": DS9, "kind": "dynamic_n1", "fault_start": 0.105}, "not a multiple"),
    ({"case": DS9, "kind": "dynamic_n1", "branches": ["1-9"]}, "no branch 1-9"),
])
def test_run_study_rejections_start_nothing(harmless, kwargs, msg):
    with pytest.raises(ValueError, match=msg):
        harmless.run_study(**kwargs)
    assert harmless.runner.jobs() == []
    assert not any(c["case"].startswith("derived/") for c in harmless.list_cases())


def test_run_study_respects_the_cpu_budget(tmp_path, monkeypatch):
    import gridpack.mcp.tools as mod
    monkeypatch.setattr(mod.os, "cpu_count", lambda: 4)
    noop = tmp_path / "noop.py"
    noop.write_text("")
    t = Tools(JobRunner(tmp_path / "jobs", max_parallel=3,
                        cli=(sys.executable, str(noop)),
                        mpiexec=(sys.executable, str(noop))))
    with pytest.raises(ValueError, match="3 jobs at a time on 2 ranks each exceed the 4 cpus"):
        t.run_study(DS9, "dynamic_n1", ranks=2)
    assert t.runner.jobs() == []


@pytest.mark.integration
def test_studies_match_runs_by_hand(tmp_path):
    t = Tools(JobRunner(tmp_path / "jobs", max_parallel=3))
    st = t.run_study(DS9, "dynamic_n1", end_time=0.3, wait_seconds=10)
    for _ in range(60):
        if not t.study_status(st["id"])["pending"]:
            break
        time.sleep(1)
    res = t.study_result(st["id"])
    assert (res["finished"], res["pending"], res["failed"]) == (9, 0, [])
    rows = res["table"]
    assert len(rows) == 9 and rows[0]["angle_spread_peak_deg"] > 1.0
    keys = [(not r["unstable"], -r["angle_spread_peak_deg"]) for r in rows]
    assert keys == sorted(keys)                           # most severe first
    # Islanding generator 3 (its step-up branch 3-6) is among the worst
    # outages even in the first 0.3 s.
    assert "3-6 BL" in {r["branch"] for r in rows[:3]}
    # One outage, redone by hand with the same edits, gives the same rows.
    item = next(r for r in rows if r["branch"] == "6-7 BL")
    hand = t.modify_case(DS9, [
        {"op": "add_fault", "from": 6, "to": 7, "ckt": "BL", "start": 0.1,
         "duration": 0.05},
        {"op": "trip_branch", "from": 6, "to": 7, "ckt": "BL", "time": 0.15}],
        replace_events=True)
    w_hand = t.job_result(t.run_dynamic_sim(hand["case"], end_time=0.3)["id"])["watch"]
    w_study = t.job_result(item["job_id"])["watch"]
    # The study watches every machine; the base case only generator 1.
    assert w_hand["columns"] == w_study["columns"][:9]
    assert [r[:9] for r in w_study["rows"]] == w_hand["rows"]

    sw = t.run_study(PF14, "load_sweep", factors=[1.0, 1.2], wait_seconds=10)
    table = t.study_result(sw["id"])["table"]
    assert [r["factor"] for r in table] == [1.0, 1.2]
    base = t.job_result(t.run_powerflow(PF14)["id"])
    assert table[0]["min_voltage"] == min(b["voltage"] for b in base["buses"])
    assert table[1]["load_mw"] == pytest.approx(1.2 * table[0]["load_mw"])
    assert table[1]["min_voltage"] < table[0]["min_voltage"]


# -------------------------------------------------------------
# observations (what HADREC records)
# -------------------------------------------------------------

def _observations(xml):
    root = ET.parse(xml).getroot()
    return [(o.findtext("type").strip(), o.findtext("busID").strip(),
             (o.findtext("generatorID") or "").strip())
            for o in root.find("Dynamic_simulation/observations")]


def test_describe_case_reports_what_a_dynamic_case_records(tools):
    assert tools.describe_case(DS9)["records"] == {
        "watched_generators": 1, "observations": {}}
    kundur = tools.describe_case("input/ds/input_kundur_two_area_wsieg1.xml")
    assert kundur["records"] == {"watched_generators": 1,
                                 "observations": {"generator": 2, "bus": 2}}
    assert "records" not in tools.describe_case(PF14)


def test_created_cases_observe_everything(tmp_path):
    d = _mine(tmp_path, "raw/9b3g.raw", "dyr/9b3g.dyr")
    t = Tools(JobRunner(tmp_path / "jobs"), cases=[d])
    c = t.create_case("mine/dyn/9b3g.raw", dyr="mine/dyn/9b3g.dyr")
    assert c["records"] == {"watched_generators": 3,
                            "observations": {"generator": 3, "bus": 9}}
    obs = _observations(d / "dyn/9b3g.xml")
    assert obs[:3] == [("generator", "1", "1"), ("generator", "2", "1"),
                       ("generator", "3", "1")]
    assert obs[3] == ("bus", "1", "") and len(obs) == 12


def test_add_observations_edit(tmp_path):
    t = Tools(JobRunner(tmp_path / "jobs"))
    d = t.modify_case(DS9, [
        {"op": "add_fault", "from": 6, "to": 7, "start": 0.1, "duration": 0.05},
        {"op": "add_observations", "generators": "all", "buses": [5, 7],
         "bus_frequency": True}], replace_events=True)
    assert d["records"] == {"watched_generators": 1, "observations": {
        "generator": 3, "bus": 2, "busfrequency": 2}}
    assert d["edits"][1]["observations"] == {"generator": 3, "bus": 2,
                                             "busfrequency": 2}
    assert d["apps"] == ["powerflow", "dsf", "hadrec"]
    xml = tmp_path / "jobs/derived" / d["case"].split("/", 1)[1]
    assert _observations(xml) == [
        ("generator", "1", "1"), ("generator", "2", "1"), ("generator", "3", "1"),
        ("bus", "5", ""), ("bus", "7", ""),
        ("busfrequency", "5", ""), ("busfrequency", "7", "")]
    # Replaces existing observations; generators by bus.
    k = t.modify_case("input/ds/input_kundur_two_area_wsieg1.xml", [
        {"op": "add_observations", "generators": [3, 4]}])
    assert k["records"]["observations"] == {"generator": 2}
    for edits, msg in [
            ([{"op": "add_observations"}], "needs generators and/or buses"),
            ([{"op": "add_observations", "buses": [50]}], "no bus 50"),
            ([{"op": "add_observations", "generators": [5]}],
             "no generator at bus 5"),
            ([{"op": "add_observations", "generators": "some"}],
             "\"all\" or a list of bus numbers"),
            ([{"op": "add_observations", "buses": [1], "bus_frequency": 1}],
             "true or false"),
            ([{"op": "add_observations", "buses": [1], "loads": "all"}],
             "unexpected loads"),
            ([{"op": "add_observations", "buses": [1]},
              {"op": "add_observations", "buses": [2]}], "one add_observations"),
            ([{"op": "add_observations", "buses": "all"}], None)]:
        if msg is None:
            with pytest.raises(ValueError, match="no Dynamic_simulation block"):
                t.modify_case(PF14, edits)
        else:
            with pytest.raises(ValueError, match=msg):
                t.modify_case(DS9, edits)
    assert len(list((tmp_path / "jobs/derived").iterdir())) == 2


def test_hadrec_refuses_a_case_that_records_nothing(harmless):
    with pytest.raises(ValueError, match="no observations.*add_observations"):
        harmless.run_hadrec(DS9)
    assert harmless.runner.jobs() == []


@pytest.mark.integration
def test_hadrec_records_the_added_observations(tmp_path):
    t = Tools(JobRunner(tmp_path / "jobs"))
    d = t.modify_case(DS9, [
        {"op": "add_observations", "generators": "all", "buses": "all"}])
    st = t.run_hadrec(d["case"], end_time=0.2)
    obs = t.job_result(st["id"])["observations"]
    # 3 generators x (speed, angle, P, Q) + 9 buses x (vmag, vangle)
    assert (st["state"], st["channels"], len(obs["columns"])) == ("done", 30, 31)
    assert obs["columns"][1:4] == ["gen_1_1_rspeed", "gen_2_1_rspeed",
                                   "gen_3_1_rspeed"]
    assert t.run_dynamic_sim(DS9, end_time=0.1)["channels"] == 8
    text = t.job_report(st["id"])
    assert "30 channels" in text and "Peak speed deviation" in text


# -------------------------------------------------------------
# create_measurements
# -------------------------------------------------------------

def test_create_measurements_needs_a_converged_power_flow(tmp_path):
    noop = tmp_path / "noop.py"
    noop.write_text("")
    t = Tools(JobRunner(tmp_path / "jobs", cli=(sys.executable, str(noop))))
    st = t.run_powerflow(PF14)              # fails: no result file
    with pytest.raises(ValueError, match="no results"):
        t.create_measurements(st["id"])
    with pytest.raises(KeyError):
        t.create_measurements("nope")
    assert not (tmp_path / "jobs/derived").exists()


@pytest.mark.integration
def test_measurements_round_trip_through_the_estimator(tmp_path):
    t = Tools(JobRunner(tmp_path / "jobs"))
    pf = t.run_powerflow(PF14)
    assert pf["case"] == PF14
    with pytest.raises(ValueError, match="is a dsf run"):
        t.create_measurements(t.run_dynamic_sim(DS9, end_time=0.05)["id"])
    # Exact meters: the estimate reproduces the power flow.
    exact = t.create_measurements(pf["id"], noise_scale=0.0, name="exact")
    assert exact["case"] == "derived/exact/exact.xml" and exact["apps"] == ["se"]
    assert exact["by_type"] == {"VM": 14, "PI": 14, "QI": 14, "PIJ": 20, "QIJ": 20}
    assert exact["files"] == {"State_estimation.networkConfiguration": "IEEE14.raw",
                              "State_estimation.measurementList": "measurements.xml"}
    se = t.run_state_estimation(exact["case"])
    assert se["converged"] is True
    cmp = t.compare_jobs(pf["id"], se["id"])
    assert abs(cmp["voltage_pu"]["max_change"]) < 2e-4 and abs(cmp["angle_deg"]["max_change"]) < 0.02
    text = t.job_report(se["id"])
    assert "Against the power flow truth: RMS voltage error 0.0000" in text
    assert "Bad data" not in text
    # Noisy meters with two planted gross errors on bus meters.
    noisy = t.create_measurements(pf["id"], seed=5, bad_data=[
        {"type": "VM", "bus": 9, "error": 0.15},
        {"type": "PI", "bus": 5, "error": 0.5}])
    assert [p["type"] for p in noisy["planted"]] == ["VM", "PI"]
    se2 = t.run_state_estimation(noisy["case"])
    text = t.job_report(se2["id"])
    assert "Bad data: 2 planted" in text and se2["converged"] is True
    truth = json.loads((tmp_path / "jobs/derived/input_14-se-1/measurements.json").read_text())
    assert truth["truth_job"] == pf["id"] and len(truth["truth_buses"]) == 14
    for kwargs, msg in [({"types": ["XX"]}, "unknown measurement type"),
                        ({"random_bad": -1}, "non-negative int"),
                        ({"noise": {"VM": -1}}, "non-negative standard"),
                        ({"noise_scale": -2}, "non-negative number"),
                        ({"name": "exact"}, "exists")]:
        with pytest.raises(ValueError, match=msg):
            t.create_measurements(pf["id"], **kwargs)


# -------------------------------------------------------------
# create_dyr and set_dyr_param
# -------------------------------------------------------------

def test_create_dyr_then_case(tmp_path):
    d = _mine(tmp_path, "raw/IEEE14.raw")
    t = Tools(JobRunner(tmp_path / "jobs"), cases=[d])
    r = t.create_dyr("mine/IEEE14.raw")
    assert r == {"dyr": "mine/IEEE14_classical.dyr", "generators": 5,
                 "models": ["GENCLS"],
                 "note": "generic textbook parameters, not validated for these "
                         "machines"}
    det = t.create_dyr("mine/IEEE14.raw", model_set="detailed", name="rich")
    assert det["dyr"] == "mine/rich.dyr" and det["models"] == ["GENROU", "SEXS", "TGOV1"]
    c = t.create_case("mine/IEEE14.raw", dyr=det["dyr"], name="dyn14")
    assert c["apps"] == ["powerflow", "dsf"] and c["records"]["watched_generators"] == 5
    for kwargs, msg in [({"raw": "mine/IEEE14.raw"}, "exists"),
                        ({"raw": "mine/IEEE14.raw", "model_set": "x"}, "model_set is one of"),
                        ({"raw": "mine/nope.raw"}, "not found"),
                        ({"raw": "mine/IEEE14.raw", "name": "a b"}, "name must be")]:
        with pytest.raises(ValueError, match=msg):
            t.create_dyr(**kwargs)
    with pytest.raises(ValueError, match="--cases DIR"):
        Tools(JobRunner(tmp_path / "jobs2")).create_dyr("x.raw")


def test_set_dyr_param_edit(tmp_path):
    t = Tools(JobRunner(tmp_path / "jobs"))
    d = t.modify_case(DS9, [
        {"op": "set_dyr_param", "bus": 2, "model": "GENCLS", "field": "H", "value": 12.5},
        {"op": "scale_load", "bus": 5, "factor": 1.1}])
    assert d["edits"][0] == {"op": "set_dyr_param", "bus": 2, "model": "GENCLS",
                             "field": "H", "value": 12.5, "was": 8.0, "now": 12.5}
    assert d["edits"][1]["now"] == {"PL": 99.0, "QL": 33.0}
    assert d["files"]["Dynamic_simulation.generatorParameters"] == "9b3g.dyr"
    from gridpack.dyredit import dyr_records
    dyr = tmp_path / "jobs/derived/input_9b3g-1/9b3g.dyr"
    assert dyr_records(dyr)[1].get("H") == 12.5 and dyr_records(dyr)[0].get("H") == 8.0
    for edits, msg in [
            ([{"op": "set_dyr_param", "bus": 9, "model": "GENCLS", "field": "H",
               "value": 1}], "no GENCLS record"),
            ([{"op": "set_dyr_param", "bus": 1, "model": "GENCLS", "field": "XD",
               "value": 1}], "has no field XD")]:
        with pytest.raises(ValueError, match=msg):
            t.modify_case(DS9, edits)
    with pytest.raises(ValueError, match="no Dynamic_simulation block"):
        t.modify_case(PF14, [{"op": "set_dyr_param", "bus": 1, "model": "GENCLS",
                              "field": "H", "value": 1}])
    assert len(list((tmp_path / "jobs/derived").iterdir())) == 1


@pytest.mark.integration
def test_synthesized_and_edited_models_run(tmp_path):
    d = _mine(tmp_path, "raw/9b3g.raw")
    t = Tools(JobRunner(tmp_path / "jobs"), cases=[d])
    # Generic detailed models on the 9-bus network, faulted, run and stable.
    dyr = t.create_dyr("mine/dyn/9b3g.raw", model_set="detailed")["dyr"]
    case = t.create_case("mine/dyn/9b3g.raw", dyr=dyr)["case"]
    f = t.modify_case(case, [{"op": "add_fault", "from": 6, "to": 7, "start": 0.1,
                              "duration": 0.05}])
    st = t.run_dynamic_sim(f["case"], end_time=1.0)
    assert (st["state"], st["exit_code"], st["channels"]) == ("done", 0, 24)
    text = t.job_report(st["id"])
    assert "Verdict:" in text and "loss of synchronism" not in text
    # A heavier machine swings less.
    base = t.run_dynamic_sim(DS9, end_time=0.5)["id"]
    heavy = t.modify_case(DS9, [{"op": "set_dyr_param", "bus": 1, "model": "GENCLS",
                                 "field": "H", "value": 16.0}])
    hv = t.run_dynamic_sim(heavy["case"], end_time=0.5)["id"]
    s0 = post_stats(t, base)
    s1 = post_stats(t, hv)
    assert s1["peak_speed_deviation"] < s0["peak_speed_deviation"]
    assert t.compare_jobs(base, hv)["first_difference_at"] <= 0.04


def post_stats(t, job_id):
    from gridpack.mcp import post
    return post.dynamic_stats(t.job_result(job_id), "dsf")
