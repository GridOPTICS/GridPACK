# -------------------------------------------------------------
# file: gridpack/mcp/tools.py
# -------------------------------------------------------------
# MCP tool functions over the job runner.
# -------------------------------------------------------------

"""Tools: named cases in, job ids and parsed ``--json`` results out.

Cases are the bundled ``data_sets`` inputs, named by their path there, and
the XMLs under each ``cases`` directory, named ``<dir name>/<path>``.  Only
names :meth:`list_cases` returns are accepted.
"""

import os
import re
import sys
import threading
from pathlib import Path
from typing import Any, Dict, List, Optional

from ..case import _OVERRIDES, Case, data_sets
from .runner import JobRunner

_APPS = ("powerflow", "ca", "dsf", "hadrec", "se", "emt")
_CASE_DIRS = ("input", "emt")
_MAX_WAIT = 10.0
_MAX_TAIL = 1000
_TIMEOUT = 600.0


class Tools:
    """The tool set; one :class:`JobRunner` holds the jobs."""

    def __init__(self, runner: Optional[JobRunner] = None, root=None,
                 timeout: float = _TIMEOUT, cases=()):
        if root is None:
            try:
                root = data_sets()
            except FileNotFoundError:
                pass
        self._root = Path(root) if root else None
        self._dirs = {}
        for d in cases:
            d = Path(d).expanduser().resolve()
            if not d.is_dir():
                raise ValueError("cases directory %s does not exist" % d)
            if d.name in _CASE_DIRS or d.name in self._dirs:
                raise ValueError("cases directory name %r is taken; rename "
                                 "the directory" % d.name)
            self._dirs[d.name] = d
        if self._root is None and not self._dirs:
            raise ValueError("no cases: no GridPACK data_sets directory "
                             "(set GRIDPACK_TEST_DATA to a GridPACK source "
                             "tree) and no cases directory")
        self.runner = runner or JobRunner()
        self.timeout = timeout
        self._cases = None
        self._lock = threading.Lock()

    # ---- cases -------------------------------------------------------

    def list_cases(self, app: Optional[str] = None) -> List[Dict[str, Any]]:
        """Available cases: applications, bus count and rank limit.

        app: only cases for this application (powerflow, ca, dsf, hadrec,
        se, emt).
        """
        cases = self._catalog()
        if app is not None:
            if app not in _APPS:
                raise ValueError("unknown application %r; expected one of %s"
                                 % (app, ", ".join(_APPS)))
            cases = [c for c in cases if app in c["apps"]]
        return [{k: v for k, v in c.items() if not k.startswith("_")}
                for c in cases]

    def describe_case(self, case: str) -> Dict[str, Any]:
        """A case's applications, data files and current solver settings."""
        entry = self._entry(case)
        c = self._load(case)
        settings = {app: {k: c.get(paths[0]) for k, paths in
                          _OVERRIDES[app].items()}
                    for app in entry["apps"]}
        return {"case": case, "apps": entry["apps"],
                "buses": entry["buses"], "max_ranks": entry["max_ranks"],
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

    def run_hadrec(self, case: str, time_step: Optional[float] = None,
                   end_time: Optional[float] = None, ranks: int = 1,
                   wait_seconds: float = _MAX_WAIT) -> Dict[str, Any]:
        """Start a HADREC dynamic run with its observations recorded; return
        its job status.  Times in seconds; waits up to wait_seconds (at
        most 10) first."""
        return self._run("hadrec", case, ranks, wait_seconds,
                         time_step=time_step, end_time=end_time)

    def run_emt(self, case: str, time_step: Optional[float] = None,
                end_time: Optional[float] = None, ranks: int = 1,
                wait_seconds: float = _MAX_WAIT) -> Dict[str, Any]:
        """Start an electromagnetic transient run; return its job status.
        Results are status only.  Times in seconds; waits up to
        wait_seconds (at most 10) first."""
        return self._run("emt", case, ranks, wait_seconds,
                         time_step=time_step, end_time=end_time)

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
        entry = self._entry(case)
        if app not in entry["apps"]:
            raise ValueError("%s cannot run %s; it supports %s"
                             % (case, app, ", ".join(entry["apps"])))
        if isinstance(ranks, int) and ranks > entry["max_ranks"]:
            raise ValueError("%s (%d buses) runs on at most %d ranks"
                             % (case, entry["buses"], entry["max_ranks"]))
        c = self._load(case)
        settings = {k: v for k, v in settings.items() if v is not None}
        if settings:
            c.override(app, **settings)
        # EMT aborts without its PETSc options file in the working dir.
        files = [entry["_petscrc"]] if app == "emt" and entry["_petscrc"] \
            else []
        job = self.runner.submit(app, c, ranks=ranks, files=files,
                                 timeout=self.timeout)
        self.runner.wait(job.id, max(0.0, min(float(wait_seconds),
                                              _MAX_WAIT)))
        return self.job_status(job.id)

    def _catalog(self):
        with self._lock:
            if self._cases is None:
                bundled = self._root / "emt/.petscrc" if self._root else None
                cases = []
                if self._root:
                    cases += _scan(self._root, _CASE_DIRS, None, bundled)
                for label, d in self._dirs.items():
                    cases += _scan(d, (".",), label, bundled)
                self._cases = cases
            return self._cases

    def _entry(self, case):
        for c in self._catalog():
            if c["case"] == case:
                return c
        raise ValueError("unknown case %r; see list_cases" % case)

    def _load(self, case):
        entry = self._entry(case)
        return _load(entry["_xml"], entry["_base"])


def _apps(case):
    """The applications a case can run, judged by its blocks."""
    root = case._root
    blocks = {el.tag for el in root}
    pf = "Powerflow" in blocks
    can = {
        "powerflow": pf,
        "ca": pf and "Contingency_analysis" in blocks,
        "dsf": pf and "Dynamic_simulation" in blocks,
        # HADREC reads the first fault event unconditionally.
        "hadrec": pf and root.find("Dynamic_simulation/Events/faultEvent")
        is not None,
        # se reads its network from its own block.
        "se": "State_estimation" in blocks,
        "emt": pf and "EMT" in blocks,
    }
    return [a for a in _APPS if can[a]]


def _load(xml, base):
    # Files are found beside the XML, then in base/raw, base/dyr, ...
    return Case.from_data_sets(xml.relative_to(base), root=base)


def _scan(base, subdirs, label, petscrc):
    """XMLs under ``base/<subdirs>`` that some application can run and whose
    data files are all present.  ``label`` marks a cases directory: its
    names get the prefix and its skips are reported on stderr."""
    out = []
    for xml in sorted(x for d in subdirs for x in (base / d).rglob("*.xml")):
        rel = xml.relative_to(base).as_posix()
        name = "%s/%s" % (label, rel) if label else rel
        # Contingency and measurement lists sit beside the inputs.
        if "<Configuration" not in xml.read_text(errors="replace"):
            continue
        try:
            c = _load(xml, base)
            files = c.files()
            apps = _apps(c)
            missing = [v for v in files.values() if not Path(v).is_absolute()]
            reason = ("data files not found: %s" % ", ".join(missing)
                      if missing else None if apps else "no application block")
        except (ValueError, OSError) as exc:
            reason = str(exc)
        if reason:
            if label:
                sys.stderr.write("gridpack-mcp: skipping %s: %s\n"
                                 % (name, reason))
            continue
        net = next((v for p, v in files.items()
                    if ".networkConfiguration" in p), None)
        buses = _count_buses(net) if net else 0
        local = xml.parent / ".petscrc"
        out.append({"case": name, "apps": apps,
                    "network": Path(net).name if net else None,
                    "buses": buses, "max_ranks": _max_ranks(buses),
                    "_xml": xml, "_base": base,
                    "_petscrc": local if local.is_file() else petscrc})
    return out


def _count_buses(raw):
    """Records in a PSS/E RAW bus section (v23 to v36)."""
    lines = Path(raw).read_text(errors="replace").splitlines()
    # v35+ put system-wide data first; older versions have 3 header lines.
    start = next((i + 1 for i, l in enumerate(lines)
                  if "BEGIN BUS DATA" in l.upper()), 3)
    n = 0
    for line in lines[start:]:
        s = line.strip()
        if not s or s.startswith("@!"):
            continue
        if re.match(r"(0|Q)(\s|/|,|$)", s):
            break
        n += 1
    return n


def _max_ranks(buses):
    # Partitioned finer, a rank gets no buses and every rank spins.
    if buses < 9:
        return 1
    if buses == 9:
        return 3
    return max(1, min(os.cpu_count() or 1, buses))
