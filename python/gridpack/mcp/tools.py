# -------------------------------------------------------------
# file: gridpack/mcp/tools.py
# -------------------------------------------------------------
# MCP tool functions over the job runner.
# -------------------------------------------------------------

"""Tools: bundled cases in, job ids and parsed ``--json`` results out.

Cases are named by their path under ``data_sets``, as :meth:`list_cases`
returns them; nothing else is accepted.
"""

import threading
from pathlib import Path
from typing import Any, Dict, List, Optional

from ..case import _OVERRIDES, Case, data_sets
from .runner import JobRunner

# Top-level block -> the application it enables.
_BLOCK_APPS = (("Powerflow", "powerflow"), ("Contingency_analysis", "ca"),
               ("Dynamic_simulation", "dsf"), ("State_estimation", "se"))
_APPS = [a for _, a in _BLOCK_APPS]
_MAX_WAIT = 10.0
_MAX_TAIL = 1000


class Tools:
    """The tool set; one :class:`JobRunner` holds the jobs."""

    def __init__(self, runner: Optional[JobRunner] = None, root=None):
        self.runner = runner or JobRunner()
        self._root = root
        self._cases = None
        self._lock = threading.Lock()

    # ---- cases -------------------------------------------------------

    def list_cases(self, app: Optional[str] = None) -> List[Dict[str, Any]]:
        """Bundled input cases and the applications each can run.

        app: only cases for this application (powerflow, ca, dsf, se).
        """
        cases = self._catalog()
        if app is not None:
            if app not in _APPS:
                raise ValueError("unknown application %r; expected one of %s"
                                 % (app, ", ".join(_APPS)))
            cases = [c for c in cases if app in c["apps"]]
        return [dict(c) for c in cases]

    def describe_case(self, case: str) -> Dict[str, Any]:
        """A case's applications, data files and current solver settings."""
        entry = self._entry(case)
        c = self._load(case)
        settings = {app: {k: c.get(paths[0]) for k, paths in
                          _OVERRIDES[app].items()}
                    for app in entry["apps"]}
        return {"case": case, "apps": entry["apps"],
                "files": {p: Path(v).name for p, v in c.files().items()},
                "settings": settings}

    # ---- runs ----------------------------------------------------------

    def run_powerflow(self, case: str, tolerance: Optional[float] = None,
                      max_iterations: Optional[int] = None, ranks: int = 1,
                      wait_seconds: float = _MAX_WAIT) -> Dict[str, Any]:
        """Start an AC power flow; return its job status.

        Waits up to wait_seconds (at most 10) for it to finish first.
        """
        return self._run("powerflow", case, ranks, wait_seconds,
                         tolerance=tolerance, max_iterations=max_iterations)

    def run_contingency(self, case: str, tolerance: Optional[float] = None,
                        max_iterations: Optional[int] = None, ranks: int = 1,
                        wait_seconds: float = _MAX_WAIT) -> Dict[str, Any]:
        """Start an N-1 contingency screen; return its job status.

        Waits up to wait_seconds (at most 10) for it to finish first.
        """
        return self._run("ca", case, ranks, wait_seconds,
                         tolerance=tolerance, max_iterations=max_iterations)

    def run_dynamic_sim(self, case: str, time_step: Optional[float] = None,
                        end_time: Optional[float] = None, ranks: int = 1,
                        wait_seconds: float = _MAX_WAIT) -> Dict[str, Any]:
        """Start a dynamic simulation; return its job status.

        time_step and end_time are in seconds.  Waits up to wait_seconds
        (at most 10) for it to finish first.
        """
        return self._run("dsf", case, ranks, wait_seconds,
                         time_step=time_step, end_time=end_time)

    def run_state_estimation(self, case: str,
                             tolerance: Optional[float] = None,
                             max_iterations: Optional[int] = None,
                             ranks: int = 1,
                             wait_seconds: float = _MAX_WAIT
                             ) -> Dict[str, Any]:
        """Start a state estimation with the case's meters; return its
        job status.  Waits up to wait_seconds (at most 10) first."""
        return self._run("se", case, ranks, wait_seconds,
                         tolerance=tolerance, max_iterations=max_iterations)

    # ---- jobs ----------------------------------------------------------

    def job_status(self, job_id: str) -> Dict[str, Any]:
        """State (queued, running, done, failed), exit code and error."""
        job = self.runner.get(job_id)
        status = job.status()
        status.pop("workdir")
        result = job.result if job.state == "done" else None
        if result is not None and "converged" in result:
            status["converged"] = result["converged"]
        return status

    def job_result(self, job_id: str) -> Dict[str, Any]:
        """The finished job's full results (the CLI's --json output)."""
        job = self.runner.get(job_id)
        if job.state in ("queued", "running"):
            raise ValueError("job %s is still %s" % (job_id, job.state))
        result = job.result
        if result is None:
            raise ValueError("job %s left no results: %s"
                             % (job_id, job.error))
        result.pop("config", None)
        return result

    def job_log(self, job_id: str, stream: str = "stdout",
                tail: int = 100) -> str:
        """The last tail lines (at most 1000) of stdout or stderr."""
        if stream not in ("stdout", "stderr"):
            raise ValueError("stream is stdout or stderr, not %r" % stream)
        return self.runner.get(job_id).log(
            stream, tail=max(1, min(int(tail), _MAX_TAIL)))

    # ---- internals -------------------------------------------------------

    def _run(self, app, case, ranks, wait_seconds, **settings):
        apps = self._entry(case)["apps"]
        if app not in apps:
            raise ValueError("%s cannot run %s; it supports %s"
                             % (case, app, ", ".join(apps)))
        c = self._load(case)
        settings = {k: v for k, v in settings.items() if v is not None}
        if settings:
            c.override(app, **settings)
        job = self.runner.submit(app, c, ranks=ranks)
        self.runner.wait(job.id, max(0.0, min(float(wait_seconds),
                                              _MAX_WAIT)))
        return self.job_status(job.id)

    def _catalog(self):
        with self._lock:
            if self._cases is None:
                self._cases = _scan(self._data())
            return self._cases

    def _entry(self, case):
        for c in self._catalog():
            if c["case"] == case:
                return c
        raise ValueError("unknown case %r; see list_cases" % case)

    def _load(self, case):
        return Case.from_data_sets(case, root=self._data())

    def _data(self):
        return self._root or data_sets()


def _scan(root):
    """Input XMLs under ``root/input`` that some application can run and
    whose data files are all present."""
    out = []
    for xml in sorted((root / "input").rglob("*.xml")):
        rel = xml.relative_to(root).as_posix()
        try:
            c = Case.from_data_sets(rel, root=root)
        except (ValueError, OSError):
            continue
        files = c.files()
        blocks = {el.tag for el in c._root}
        # se reads its network from its own block, the rest from Powerflow.
        apps = [a for b, a in _BLOCK_APPS
                if b in blocks and (a == "se" or "Powerflow" in blocks)]
        if not apps or not all(Path(v).is_absolute() for v in files.values()):
            continue
        net = next((v for p, v in files.items()
                    if ".networkConfiguration" in p), None)
        out.append({"case": rel, "apps": apps,
                    "network": Path(net).name if net else None})
    return out
