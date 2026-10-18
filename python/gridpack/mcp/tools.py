# -------------------------------------------------------------
# file: gridpack/mcp/tools.py
# -------------------------------------------------------------
# MCP tool functions over the job runner.
# -------------------------------------------------------------

"""Tools: named cases in, job ids and parsed ``--json`` results out.

Cases are the bundled ``data_sets`` inputs, named by their path there, and
the XMLs under each ``cases`` directory, named ``<dir name>/<path>``.  Only
names :meth:`list_cases` returns are accepted; :meth:`create_case` adds one
for a RAW/DYR pair under a cases directory.
"""

import json
import os
import re
import shutil
import sys
import threading
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any, Dict, List, Optional

from ..case import _OVERRIDES, Case, data_sets
from ..rawedit import OPS, edit_raw
from .events import EVENT_OPS, append_events, plan_events
from . import plots, post
from .runner import JobRunner
from .templates import build_xml, raw_generators, raw_version

_APPS = ("powerflow", "ca", "dsf", "hadrec", "se", "emt")
_CASE_DIRS = ("input", "emt")
_DERIVED = "derived"
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
            if d.name in _CASE_DIRS or d.name == _DERIVED or d.name in self._dirs:
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

    def create_case(self, raw: str, dyr: Optional[str] = None,
                    name: Optional[str] = None) -> Dict[str, Any]:
        """Write an input XML for a RAW (and DYR) file beside it under a
        cases directory and list it; return its description.

        raw, dyr: "<dir name>/<path>" as list_cases names things, or a
        path under a cases directory.  name: the XML's stem (default:
        the RAW's).  The case runs powerflow, and dsf with a DYR (every
        generator watched, no events yet).
        """
        if not self._dirs:
            raise ValueError("create_case needs a cases directory: start "
                             "gridpack-mcp with --cases DIR")
        raw_path = self._under_cases(raw, "RAW")
        dyr_path = self._under_cases(dyr, "DYR") if dyr else None
        stem = name or raw_path.stem
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", stem):
            raise ValueError("name must be letters, digits, _ . or -, "
                             "not %r" % stem)
        xml = raw_path.parent / (stem + ".xml")
        if xml.exists():
            raise ValueError("%s exists; pick another name" % self._name(xml))
        gens = raw_generators(raw_path) if dyr_path else []
        if dyr_path and not gens:
            raise ValueError("no generator records in %s" % raw_path.name)
        xml.write_text(build_xml(
            raw_path.name, raw_version(raw_path),
            dyr=os.path.relpath(dyr_path, xml.parent) if dyr_path else None,
            generators=gens))
        with self._lock:
            self._cases = None
        try:
            return self.describe_case(self._name(xml))
        except ValueError:
            xml.unlink()
            raise

    def modify_case(self, case: str, edits: List[Dict[str, Any]],
                    name: Optional[str] = None) -> Dict[str, Any]:
        """Derive a case from ``case`` with its network edited and events
        added; list it as derived/<name>/<name>.xml and return its
        description plus what each edit changed.

        edits: a list of {"op": ..., ...}; ids default to the only one
        there is.  Network ops: set_load_active_power {bus, mw};
        set_load_reactive_power {bus, mvar}; scale_load {bus, factor};
        set_load_in_service {bus, in_service}; set_generator_active_power
        {bus, mw}; set_generator_voltage_setpoint {bus, vm_pu};
        set_generator_in_service {bus, in_service}; set_shunt {bus, mvar,
        mw}; set_shunt_in_service {bus, in_service}; set_branch_in_service
        {from, to, ckt, in_service}; set_param {kind, bus|from+to+ckt, id,
        field, value} with a PSS/E RAW column name.  Event ops (seconds,
        on a time step): add_fault {from, to, start, duration} (solid
        fault at the from bus); add_bus_fault {bus, start, duration, g_pu,
        b_pu}; trip_branch {from, to, ckt, time}; trip_generator {bus, id,
        time}.  hadrec needs a fault as the first event.  Everything else
        in the case is inherited.
        """
        entry = self._entry(case)
        if not isinstance(edits, list) or not edits:
            raise ValueError("edits must be a non-empty list")
        for n, e in enumerate(edits, 1):
            if not isinstance(e, dict) or "op" not in e:
                raise ValueError("edit %d must be an object with an \"op\"" % n)
            if e["op"] not in OPS and e["op"] not in EVENT_OPS:
                raise ValueError("edit %d: unknown op %r; one of %s" % (
                    n, e["op"], ", ".join(list(OPS) + list(EVENT_OPS))))
        static = [e for e in edits if e["op"] in OPS]
        dynamic = [e for e in edits if e["op"] in EVENT_OPS]
        c = self._load(case)
        tag = next((p for p in c.files() if ".networkConfiguration" in p), None)
        if tag is None:
            raise ValueError("%s has no network file to edit" % case)
        src = Path(c.files()[tag])
        events = []
        if dynamic:
            if "dsf" not in entry["apps"]:
                raise ValueError("%s has no Dynamic_simulation block to hold "
                                 "events" % case)
            events = plan_events(dynamic, src,
                                 float(c.get("Dynamic_simulation.timeStep")),
                                 float(c.get("Dynamic_simulation.simulationTime")))
        root = self.runner.root / _DERIVED
        stem = name or self._free_name(root, Path(entry["_xml"]).stem)
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", stem):
            raise ValueError("name must be letters, digits, _ . or -, "
                             "not %r" % stem)
        d = root / stem
        if d.exists():
            raise ValueError("%s/%s exists; pick another name"
                             % (_DERIVED, stem))
        d.mkdir(parents=True)
        try:
            changed = iter(edit_raw(src, d / src.name, static) if static else [])
            if static:
                c.set(tag, str(d / src.name))
            xml = c.write(d, name=stem + ".xml")
            added = {id(e): ET.tostring(el, encoding="unicode")
                     for e, el in events}
            if events:
                append_events(xml, [el for _, el in events])
            report = [dict(e, xml=added[id(e)]) if id(e) in added
                      else next(changed) for e in edits]
            (d / "edits.json").write_text(json.dumps(
                {"base": case, "edits": report}, indent=1))
        except (ValueError, OSError):
            shutil.rmtree(d, ignore_errors=True)
            raise
        with self._lock:
            self._dirs.setdefault(_DERIVED, root)
            self._cases = None
        out = self.describe_case("%s/%s/%s" % (_DERIVED, stem, xml.name))
        out.update(base=case, edits=report)
        return out

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

    # ---- post-processing -------------------------------------------------

    def job_series(self, job_id: str, columns: Optional[List[str]] = None,
                   max_points: int = 200) -> Dict[str, Any]:
        """Chosen channels of a finished dynamic run (dsf watch, hadrec
        observations), decimated to at most max_points rows.

        columns: channel names or globs such as "*_speed" (time is always
        the first column); None means every channel.
        """
        job, result = self._finished(job_id)
        return post.series(result, job.app, columns,
                           max(2, min(int(max_points), 100000)))

    def job_report(self, job_id: str) -> str:
        """A short text verdict on a finished job: convergence, violations,
        totals, events, speed deviation, damping and angle spread."""
        job, result = self._finished(job_id)
        return post.report(result, job.app, job.workdir)

    def compare_jobs(self, job_a: str, job_b: str) -> Dict[str, Any]:
        """What changed from job_a to job_b: bus voltage and angle changes,
        totals and branch differences for static runs; per-channel largest
        change and the time it happened for dynamic runs."""
        ja, ra = self._finished(job_a)
        jb, rb = self._finished(job_b)
        return post.compare(ra, rb, ja.app, jb.app)

    def plot_job(self, job_id: str, kind: Optional[str] = None,
                 channels: Optional[List[str]] = None,
                 title: Optional[str] = None) -> bytes:
        """A PNG of a finished job.  Dynamic runs: kind "series" over
        channels (names or globs; default every generator speed), one
        panel per quantity, events shaded.  Power flow / state estimation:
        kind "voltage" (default) or "angle" over the buses.  Needs
        matplotlib (pip install 'gridpack[plot]')."""
        job, result = self._finished(job_id)
        events = post.events_from_xml(job.workdir) if job.app in ("dsf", "hadrec") \
            else ()
        return plots.render(result, job.app, kind, channels, title, events)

    # ---- internals -------------------------------------------------------

    def _finished(self, job_id):
        job = self.runner.get(job_id)
        if job.state in ("queued", "running"):
            raise ValueError("job %s is still %s" % (job_id, job.state))
        result = job.result
        if result is None:
            raise ValueError("job %s left no results: %s" % (job_id, job.error))
        return job, result

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

    def _under_cases(self, spec, kind):
        """The file ``spec`` names inside a cases directory, resolved."""
        p = Path(spec).expanduser()
        candidates = [p] if p.is_absolute() else []
        if p.parts and p.parts[0] in self._dirs:
            candidates.append(self._dirs[p.parts[0]].joinpath(*p.parts[1:]))
        for c in candidates:
            c = c.resolve()
            if any(_under(c, d) for d in self._dirs.values()):
                if not c.is_file():
                    raise ValueError("%s file %s not found" % (kind, spec))
                return c
        raise ValueError("%s file must be under a cases directory (%s), as "
                         "<dir name>/<path>; got %r"
                         % (kind, ", ".join(self._dirs), spec))

    @staticmethod
    def _free_name(root, base):
        n = 1
        while (root / ("%s-%d" % (base, n))).exists():
            n += 1
        return "%s-%d" % (base, n)

    def _name(self, path):
        """The catalog name of a file under a cases directory."""
        for label, d in self._dirs.items():
            if _under(path, d):
                return "%s/%s" % (label, path.relative_to(d).as_posix())
        raise ValueError("%s is not under a cases directory" % path)


def _apps(case):
    """The applications a case can run, judged by its blocks."""
    root = case._root
    blocks = {el.tag for el in root}
    pf = "Powerflow" in blocks
    can = {
        "powerflow": pf,
        "ca": pf and "Contingency_analysis" in blocks,
        "dsf": pf and "Dynamic_simulation" in blocks,
        # HADREC initialises from the first event, which must be a fault.
        "hadrec": pf and _first_event(root) in ("faultEvent", "BusFault"),
        # se reads its network from its own block.
        "se": "State_estimation" in blocks,
        "emt": pf and "EMT" in blocks,
    }
    return [a for a in _APPS if can[a]]


def _first_event(root):
    events = root.find("Dynamic_simulation/Events")
    return events[0].tag if events is not None and len(events) else None


def _under(path, directory):
    try:
        path.relative_to(directory)
        return True
    except ValueError:
        return False


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
