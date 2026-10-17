# -------------------------------------------------------------
# file: python/tests/test_runner.py
# -------------------------------------------------------------
# gridpack.mcp.runner: job states, errors, timeouts, cancel.
# -------------------------------------------------------------

import json
import os
import shutil
import sys
import time

import pytest

from gridpack import Case
from gridpack.mcp.runner import JobRunner

# Behaves as fake.json in its cwd says: sleep, spawn a grandchild, write
# stderr / result.json, exit.
_FAKE_CLI = r'''
import json, os, subprocess, sys, time
spec = json.load(open("fake.json"))
open("argv.json", "w").write(json.dumps(sys.argv[1:]))
if spec.get("grandchild"):
    p = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(60)"])
    open("grandchild.pid", "w").write(str(p.pid))
time.sleep(spec.get("sleep", 0))
sys.stderr.write(spec.get("stderr", ""))
if "result" in spec:
    open("result.json", "w").write(json.dumps(spec["result"]))
print("fake stdout")
sys.exit(spec.get("exit", 0))
'''

# Records its argv, drops "-np N" and runs the rest.
_FAKE_MPIEXEC = r'''
import subprocess, sys
open("mpiexec_argv.txt", "w").write(" ".join(sys.argv[1:3]))
sys.exit(subprocess.call(sys.argv[3:]))
'''


@pytest.fixture
def fake(tmp_path):
    (tmp_path / "cli.py").write_text(_FAKE_CLI)
    (tmp_path / "mpiexec.py").write_text(_FAKE_MPIEXEC)
    runner = JobRunner(tmp_path / "jobs",
                       cli=(sys.executable, str(tmp_path / "cli.py")),
                       mpiexec=(sys.executable, str(tmp_path / "mpiexec.py")))

    def submit(spec, **kw):
        f = tmp_path / ("fake-%d.json" % len(runner.jobs()))
        f.write_text(json.dumps(spec))
        dest = tmp_path / "fake.json"
        shutil.copyfile(f, dest)
        case = Case.from_data_sets("input/powerflow/input_14.xml")
        return runner.submit(kw.pop("app", "powerflow"), case,
                             files=[dest], **kw)

    return runner, submit


def test_done_with_result_and_logs(fake):
    runner, submit = fake
    job = runner.wait(submit({"result": {"converged": True}}).id, 60)
    assert (job.state, job.exit_code, job.error) == ("done", 0, None)
    assert job.result == {"converged": True}
    assert job.log() == "fake stdout"
    assert json.loads((job.workdir / "argv.json").read_text()) == [
        "powerflow", "input_14.xml", "--no-timer", "--json", "result.json"]
    assert sorted(p.name for p in job.workdir.glob("*.raw")) == ["IEEE14.raw"]


def test_not_converged_is_a_finished_job(fake):
    runner, submit = fake
    job = runner.wait(submit({"exit": 3, "result": {"converged": False}}).id,
                      60)
    assert (job.state, job.exit_code) == ("done", 3)


@pytest.mark.parametrize("spec, error", [
    ({"exit": 1, "result": {"error": "ValueError: bad"}}, "ValueError: bad"),
    ({"exit": 2, "result": {"error": None}, "stderr": "Error: no list\n"},
     "Error: no list"),
    ({"exit": 139, "stderr": "Segmentation fault\n"},
     "exit 139 without a result file: Segmentation fault"),
])
def test_failures_carry_the_reason(fake, spec, error):
    runner, submit = fake
    job = runner.wait(submit(spec).id, 60)
    assert (job.state, job.error) == ("failed", error)


def test_ranks_go_through_mpiexec(fake):
    runner, submit = fake
    one = runner.wait(submit({"result": {}}).id, 60)
    two = runner.wait(submit({"result": {}}, ranks=2).id, 60)
    assert not (one.workdir / "mpiexec_argv.txt").exists()
    assert (two.workdir / "mpiexec_argv.txt").read_text() == "-np 2"
    assert two.state == "done" and one.workdir != two.workdir


def test_timeout_kills_the_process_group(fake):
    runner, submit = fake
    job = runner.wait(submit({"sleep": 60, "grandchild": True},
                             timeout=2).id, 60)
    assert (job.state, job.error) == ("failed", "timed out after 2 s")
    pid = int((job.workdir / "grandchild.pid").read_text())
    time.sleep(0.5)
    with pytest.raises(ProcessLookupError):
        os.kill(pid, 0)


def test_cancel_running_and_queued(fake):
    runner, submit = fake
    first = submit({"sleep": 60})          # max_parallel=1: holds the slot
    second = submit({"result": {}})
    deadline = time.time() + 30
    while first.state != "running" and time.time() < deadline:
        time.sleep(0.05)
    assert second.state == "queued"
    runner.cancel(second.id)
    runner.cancel(first.id)
    for job in (runner.wait(first.id, 30), runner.wait(second.id, 30)):
        assert (job.state, job.error) == ("failed", "cancelled")
    # Never launched, not launched and killed.
    assert not (second.workdir / "argv.json").exists()
    assert not (second.workdir / "stdout.txt").exists()


@pytest.mark.parametrize("app, ranks", [("pf", 1), ("powerflow", 0),
                                        ("powerflow", True)])
def test_submit_validates(fake, app, ranks):
    runner, submit = fake
    with pytest.raises(ValueError):
        submit({}, app=app, ranks=ranks)
    assert runner.jobs() == [] and list(runner.root.iterdir()) == []


@pytest.mark.integration
@pytest.mark.parametrize("ranks", [1, 2])
def test_real_powerflow(tmp_path, ranks):
    if ranks > 1 and shutil.which("mpiexec") is None:
        pytest.skip("mpiexec not on PATH")
    runner = JobRunner(tmp_path)
    case = Case.from_data_sets("input/powerflow/input_14.xml")
    job = runner.wait(runner.submit("powerflow", case, ranks=ranks).id, 300)
    assert (job.state, job.exit_code) == ("done", 0), job.log("stderr")
    assert job.result["converged"] is True and len(job.result["buses"]) == 14
