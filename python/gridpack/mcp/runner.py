# -------------------------------------------------------------
# file: gridpack/mcp/runner.py
# -------------------------------------------------------------
# Run gridpack CLI jobs in subprocesses and track them.
# -------------------------------------------------------------

"""Background jobs: one ``gridpack <app> ... --json`` subprocess each.

Each job gets a private directory holding the staged case, ``stdout.txt``,
``stderr.txt`` and ``result.json``.  Jobs stay in memory only.
"""

import json
import os
import shutil
import signal
import subprocess
import tempfile
import threading
import time
import uuid
from pathlib import Path

__all__ = ["Job", "JobRunner", "APPS"]

APPS = ("powerflow", "ca", "dsf", "hadrec", "se", "emt")

# CLI exit codes; 3 is a finished analysis that did not converge.
_FINISHED = (0, 3)


class Job:
    """One run.  ``state`` is queued, running, done or failed."""

    def __init__(self, job_id, app, ranks, workdir, timeout):
        self.id = job_id
        self.app = app
        self.ranks = ranks
        self.workdir = workdir
        self.timeout = timeout
        self.state = "queued"
        self.exit_code = None
        self.error = None
        self.command = None
        self.submitted = time.time()
        self.started = None
        self.finished = None
        self._proc = None
        self._cancelled = False
        self._done = threading.Event()

    @property
    def result(self):
        """The parsed ``--json`` output, or None."""
        try:
            return json.loads((self.workdir / "result.json").read_text())
        except (OSError, ValueError):
            return None

    def log(self, stream="stdout", tail=None):
        """The last ``tail`` lines of ``stdout`` or ``stderr``."""
        path = self.workdir / ("%s.txt" % stream)
        lines = path.read_text(errors="replace").splitlines() \
            if path.is_file() else []
        return "\n".join(lines[-tail:] if tail else lines)

    def status(self):
        """A JSON-ready summary."""
        return {"id": self.id, "app": self.app, "ranks": self.ranks,
                "state": self.state, "exit_code": self.exit_code,
                "error": self.error, "workdir": str(self.workdir),
                "submitted": self.submitted, "started": self.started,
                "finished": self.finished}


class JobRunner:
    """Queue jobs and run up to ``max_parallel`` at a time.

    ``cli`` and ``mpiexec`` are the commands used, for tests and installs
    where they are not on PATH.
    """

    def __init__(self, root=None, *, max_parallel=1, cli=("gridpack",),
                 mpiexec=("mpiexec",)):
        self.root = Path(root or tempfile.mkdtemp(prefix="gridpack-mcp-"))
        self.root.mkdir(parents=True, exist_ok=True)
        self._cli = list(cli)
        self._mpiexec = list(mpiexec)
        self.max_parallel = max_parallel
        self._slots = threading.Semaphore(max_parallel)
        self._jobs = {}
        self._lock = threading.Lock()

    def submit(self, app, case, *, ranks=1, timeout=None, files=()):
        """Stage ``case`` (a :class:`gridpack.Case`) and queue a run.

        ``files`` are copied into the job directory as well.
        """
        if app not in APPS:
            raise ValueError("unknown application %r; expected one of %s"
                             % (app, ", ".join(APPS)))
        if not isinstance(ranks, int) or isinstance(ranks, bool) or ranks < 1:
            raise ValueError("ranks must be a positive int, not %r" % (ranks,))
        job_id = uuid.uuid4().hex[:12]
        workdir = self.root / job_id
        xml = case.write(workdir)
        for f in files:
            shutil.copyfile(f, workdir / Path(f).name)
        job = Job(job_id, app, ranks, workdir, timeout)
        cmd = self._cli + [app, xml.name, "--no-timer", "--json",
                           "result.json"]
        job.command = (self._mpiexec + ["-np", str(ranks)] + cmd
                       if ranks > 1 else cmd)
        with self._lock:
            self._jobs[job_id] = job
        threading.Thread(target=self._run, args=(job,), daemon=True).start()
        return job

    def get(self, job_id):
        with self._lock:
            if job_id not in self._jobs:
                raise KeyError("no job %r" % job_id)
            return self._jobs[job_id]

    def jobs(self):
        with self._lock:
            return list(self._jobs.values())

    def wait(self, job_id, timeout=None):
        """Block until the job ends or ``timeout`` passes; return it."""
        job = self.get(job_id)
        job._done.wait(timeout)
        return job

    def cancel(self, job_id):
        """Stop a queued or running job.  No-op once it has ended."""
        job = self.get(job_id)
        with self._lock:
            job._cancelled = True
            proc = job._proc
        if proc is not None and proc.poll() is None:
            _kill(proc)
        return job

    def _run(self, job):
        try:
            self._execute(job)
        except Exception as exc:          # never leave a job "running"
            self._end(job, "failed", job.exit_code,
                      "%s: %s" % (type(exc).__name__, exc))

    def _execute(self, job):
        with self._slots:
            with self._lock:
                cancelled = job._cancelled
                if not cancelled:
                    job.state, job.started = "running", time.time()
            if cancelled:
                return self._end(job, "failed", error="cancelled")
            try:
                with open(job.workdir / "stdout.txt", "w") as out, \
                        open(job.workdir / "stderr.txt", "w") as err:
                    # Own process group, so a kill reaches mpiexec's ranks.
                    proc = subprocess.Popen(
                        job.command, cwd=job.workdir, stdout=out, stderr=err,
                        stdin=subprocess.DEVNULL, start_new_session=True)
                    with self._lock:
                        job._proc = proc
                        cancelled = job._cancelled
                    if cancelled:
                        _kill(proc)
                    try:
                        rc = proc.wait(timeout=job.timeout)
                    except subprocess.TimeoutExpired:
                        _kill(proc)
                        return self._end(job, "failed", proc.wait(),
                                         "timed out after %s s" % job.timeout)
            except OSError as exc:
                return self._end(job, "failed", error="could not start: %s"
                                 % exc)
        if job._cancelled:
            return self._end(job, "failed", rc, "cancelled")
        result = job.result
        if result is None:
            return self._end(job, "failed", rc,
                             "exit %d without a result file: %s"
                             % (rc, job.log("stderr", tail=1)))
        if rc in _FINISHED:
            return self._end(job, "done", rc)
        return self._end(job, "failed", rc, result.get("error")
                         or job.log("stderr", tail=1) or "exit %d" % rc)

    def _end(self, job, state, exit_code=None, error=None):
        with self._lock:
            job.state, job.exit_code, job.error = state, exit_code, error
            job.finished = time.time()
        job._done.set()


def _kill(proc):
    try:
        os.killpg(proc.pid, signal.SIGTERM)
        proc.wait(timeout=5)
    except subprocess.TimeoutExpired:
        os.killpg(proc.pid, signal.SIGKILL)
    except ProcessLookupError:
        pass
