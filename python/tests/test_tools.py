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
