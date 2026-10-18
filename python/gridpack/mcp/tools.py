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
import time
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any, Dict, List, Optional

from ..case import _OVERRIDES, Case, data_sets
from ..dyredit import MODEL_SETS, edit_dyr, synthesize_dyr
from ..rawedit import OPS, edit_raw
from .events import (EVENT_OPS, append_events, plan_events, raw_elements,
                     set_observations, watch_all_generators)
from . import measurements, plots, post, study
from .runner import JobRunner
from .templates import build_xml, raw_generators, raw_version

_APPS = ("powerflow", "ca", "dsf", "hadrec", "se", "emt")
_CASE_DIRS = ("input", "emt")
_DERIVED = "derived"
_OBSERVE = "add_observations"
_DYR = "set_dyr_param"
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
        self._studies = {}
        self._job_case = {}
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
        out = {"case": case, "apps": entry["apps"],
               "buses": entry["buses"], "max_ranks": entry["max_ranks"],
               "files": {p: Path(v).name for p, v in c.files().items()},
               "settings": settings}
        if "dsf" in entry["apps"]:
            out["records"] = _records(c._root)
        return out

    def create_case(self, raw: str, dyr: Optional[str] = None,
                    name: Optional[str] = None) -> Dict[str, Any]:
        """Write an input XML for a RAW (and DYR) file beside it under a
        cases directory and list it; return its description.

        raw, dyr: "<dir name>/<path>" as list_cases names things, or a
        path under a cases directory.  name: the XML's stem (default:
        the RAW's).  The case runs powerflow, and dsf with a DYR (every
        generator watched and observed, every bus observed, no events yet).
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
            generators=gens,
            buses=sorted(raw_elements(raw_path)[0]) if dyr_path else ()))
        with self._lock:
            self._cases = None
        try:
            return self.describe_case(self._name(xml))
        except ValueError:
            xml.unlink()
            raise

    def create_dyr(self, raw: str, model_set: str = "classical",
                   name: Optional[str] = None) -> Dict[str, Any]:
        """Write a DYR beside a RAW under a cases directory, giving every
        generator the generic models of model_set: "classical" (GENCLS) or
        "detailed" (GENROU + SEXS + TGOV1), at textbook parameters.

        The parameters are not those of any real machine: use the file for
        teaching, screening and comparing scenarios, and tune it with
        modify_case set_dyr_param.  Pass the result to create_case as dyr.
        """
        if not self._dirs:
            raise ValueError("create_dyr needs a cases directory: start "
                             "gridpack-mcp with --cases DIR")
        raw_path = self._under_cases(raw, "RAW")
        if model_set not in MODEL_SETS:
            raise ValueError("model_set is one of %s, not %r"
                             % (", ".join(MODEL_SETS), model_set))
        stem = name or "%s_%s" % (raw_path.stem, model_set)
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", stem):
            raise ValueError("name must be letters, digits, _ . or -, "
                             "not %r" % stem)
        dyr = raw_path.parent / (stem + ".dyr")
        if dyr.exists():
            raise ValueError("%s exists; pick another name" % self._name(dyr))
        gens = raw_generators(raw_path)
        if not gens:
            raise ValueError("no generator records in %s" % raw_path.name)
        dyr.write_text(synthesize_dyr(gens, model_set))
        return {"dyr": self._name(dyr), "generators": len(gens),
                "models": [m for m, _ in MODEL_SETS[model_set]],
                "note": "generic textbook parameters, not validated for "
                        "these machines"}

    def modify_case(self, case: str, edits: List[Dict[str, Any]],
                    name: Optional[str] = None,
                    replace_events: bool = False) -> Dict[str, Any]:
        """Derive a case from ``case`` with its network edited and events
        added (replacing the base case's events if replace_events); list
        it as derived/<name>/<name>.xml and return its description plus
        what each edit changed.

        edits: a list of {"op": ..., ...}; ids default to the only one
        there is.  Network ops: set_load_active_power {bus, mw};
        set_load_reactive_power {bus, mvar}; scale_load {bus, factor};
        set_load_in_service {bus, in_service}; set_generator_active_power
        {bus, mw}; set_generator_voltage_setpoint {bus, vm_pu};
        set_generator_in_service {bus, in_service}; set_shunt {bus, mvar,
        mw}; set_shunt_in_service {bus, in_service}; set_branch_in_service
        {from, to, ckt, in_service}; scale_all_loads {factor};
        scale_all_generation {factor}; set_param {kind, bus or from+to+ckt,
        id, field, value} with a PSS/E RAW column name.  Event ops (seconds,
        on a time step): add_fault {from, to, start, duration} (solid
        fault at the from bus); add_bus_fault {bus, start, duration, g_pu,
        b_pu}; trip_branch {from, to, ckt, time}; trip_generator {bus, id,
        time}.  hadrec needs a fault as the first event and observations:
        add_observations {generators: "all" or [bus...], buses: "all" or
        [bus...], bus_frequency} sets what it records (replacing any).
        set_dyr_param {bus, id, model, field, value} changes a dynamic model
        parameter in the case's DYR (fields named as GridPACK's parsers do,
        e.g. GENROU H, XDP; SEXS K; TGOV1 R; or p3, p4, ... by position).
        Everything else in the case is inherited.
        """
        entry = self._entry(case)
        if not isinstance(edits, list) or not edits:
            raise ValueError("edits must be a non-empty list")
        for n, e in enumerate(edits, 1):
            if not isinstance(e, dict) or "op" not in e:
                raise ValueError("edit %d must be an object with an \"op\"" % n)
            if e["op"] not in OPS and e["op"] not in EVENT_OPS \
                    and e["op"] not in (_OBSERVE, _DYR):
                raise ValueError("edit %d: unknown op %r; one of %s" % (
                    n, e["op"], ", ".join(list(OPS) + list(EVENT_OPS)
                                          + [_OBSERVE, _DYR])))
        static = [e for e in edits if e["op"] in OPS]
        dynamic = [e for e in edits if e["op"] in EVENT_OPS]
        observe = [e for e in edits if e["op"] == _OBSERVE]
        dyr_edits = [e for e in edits if e["op"] == _DYR]
        if len(observe) > 1:
            raise ValueError("give one add_observations edit")
        if observe:
            _check_observe(observe[0])
        c = self._load(case)
        tag = next((p for p in c.files() if ".networkConfiguration" in p), None)
        if tag is None:
            raise ValueError("%s has no network file to edit" % case)
        src = Path(c.files()[tag])
        events = []
        if (dynamic or observe or dyr_edits) and "dsf" not in entry["apps"]:
            raise ValueError("%s has no Dynamic_simulation block to hold "
                             "events, observations or dynamic models" % case)
        dyr_tag = next((p for p in c.files() if p.endswith(".generatorParameters")),
                       None) if dyr_edits else None
        if dyr_edits and dyr_tag is None:
            raise ValueError("%s names no DYR file to edit" % case)
        if dynamic:
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
            dyr_report = iter([])
            if dyr_edits:
                dyr_src = Path(c.files()[dyr_tag])
                dyr_report = iter(edit_dyr(dyr_src, d / dyr_src.name, [
                    {k: v for k, v in e.items() if k != "op"} for e in dyr_edits]))
                c.set(dyr_tag, str(d / dyr_src.name))
            xml = c.write(d, name=stem + ".xml")
            added = {id(e): ET.tostring(el, encoding="unicode")
                     for e, el in events}
            if events or replace_events:
                append_events(xml, [el for _, el in events], replace_events)
            for e in observe:
                added[id(e)] = set_observations(
                    xml, src, e.get("generators"), e.get("buses"),
                    e.get("bus_frequency", False))
            report = [dict(e, observations=added[id(e)]) if e["op"] == _OBSERVE
                      else dict(next(dyr_report), op=_DYR) if e["op"] == _DYR
                      else dict(e, xml=added[id(e)]) if id(e) in added
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

    def create_measurements(self, job_id: str, types: Optional[List[str]] = None,
                            noise: Optional[Dict[str, float]] = None,
                            noise_scale: float = 1.0, coverage: float = 1.0,
                            bad_data: Optional[List[Dict[str, Any]]] = None,
                            random_bad: int = 0, bad_sigma: float = 20.0,
                            seed: int = 0, name: Optional[str] = None
                            ) -> Dict[str, Any]:
        """Make a state estimation case from a finished power flow job:
        meters taken from its solution with noise and optional bad data,
        listed as derived/<name>/<name>.xml; returns its description, the
        meter counts and what was planted.

        types: VM, VA (deg), PI, QI, PIJ, QIJ, PJI, QJI (p.u.); default VM,
        PI, QI, PIJ, QIJ on every bus and branch.  noise: std per type
        (defaults VM 0.004, PI/QI 0.01, flows 0.008, VA 0.1), scaled by
        noise_scale (0 = exact).  coverage: fraction of meters kept.
        bad_data: [{type, bus | from, to, ckt, error}] adds error to those
        meters; random_bad plants that many errors of bad_sigma standard
        deviations.  seed makes it repeatable.  The SE job's job_report
        then compares the estimate with this truth.
        """
        job, result = self._finished(job_id)
        if job.app != "powerflow":
            raise ValueError("job %s is a %s run; measurements come from a "
                             "power flow" % (job_id, job.app))
        if not result.get("converged"):
            raise ValueError("job %s did not converge; its solution is not a "
                             "truth" % job_id)
        base = self._job_case.get(job_id)
        if base is None:
            raise ValueError("job %s has no case on record" % job_id)
        if not isinstance(random_bad, int) or isinstance(random_bad, bool) \
                or random_bad < 0:
            raise ValueError("random_bad must be a non-negative int")
        for k, v in (("noise_scale", noise_scale), ("bad_sigma", bad_sigma)):
            if isinstance(v, bool) or not isinstance(v, (int, float)) or v < 0:
                raise ValueError("%s must be a non-negative number" % k)
        if noise is not None:
            if not isinstance(noise, dict) or any(
                    t not in measurements.TYPES or isinstance(v, bool)
                    or not isinstance(v, (int, float)) or v < 0
                    for t, v in noise.items()):
                raise ValueError("noise maps measurement types to non-negative "
                                 "standard deviations")
        c = self._load(base)
        tag = next((p for p in c.files() if ".networkConfiguration" in p), None)
        if tag is None:
            raise ValueError("%s has no network file" % base)
        raw = Path(c.files()[tag])
        meters, summary = measurements.synthesize(
            result, measurements.raw_sbase(raw), types, noise, noise_scale,
            coverage, bad_data or (), random_bad, bad_sigma, seed)
        root = self.runner.root / _DERIVED
        stem = name or self._free_name(root, Path(self._entry(base)["_xml"]).stem + "-se")
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", stem):
            raise ValueError("name must be letters, digits, _ . or -, "
                             "not %r" % stem)
        d = root / stem
        if d.exists():
            raise ValueError("%s/%s exists; pick another name" % (_DERIVED, stem))
        d.mkdir(parents=True)
        try:
            measurements.write_measurements(d / "measurements.xml", meters)
            (d / (stem + ".xml")).write_text(measurements.se_xml(
                str(raw), raw_version(raw), "measurements.xml"))
            xml = Case(d / (stem + ".xml")).write(d)   # copies the RAW beside
            measurements.dump(d / "measurements.json", result,
                              dict(summary, truth_job=job_id, base=base))
        except (ValueError, OSError):
            shutil.rmtree(d, ignore_errors=True)
            raise
        with self._lock:
            self._dirs.setdefault(_DERIVED, root)
            self._cases = None
        out = self.describe_case("%s/%s/%s" % (_DERIVED, stem, xml.name))
        out.update(truth_job=job_id, base=base, **summary)
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
        status["case"] = self._job_case.get(job_id)
        result = job.result if job.state == "done" else None
        if result is not None and "converged" in result:
            status["converged"] = result["converged"]
        t = post.table(result, job.app) if result is not None else None
        if job.app in ("dsf", "hadrec") and result is not None:
            status["channels"] = len(t[0]) - 1 if t else 0
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

    # ---- studies ---------------------------------------------------------

    def run_study(self, case: str, kind: str,
                  factors: Optional[List[float]] = None,
                  branches: Optional[List[str]] = None, ranks: int = 1,
                  fault_start: float = 0.1, fault_duration: float = 0.05,
                  end_time: Optional[float] = None,
                  scale_generation: bool = False,
                  wait_seconds: float = _MAX_WAIT) -> Dict[str, Any]:
        """Fan a study out as one derived case and one job per point;
        return its status (see study_status / study_result).

        kind "load_sweep": each factor scales every load (and every
        generator with scale_generation) and runs a power flow.  kind
        "dynamic_n1": each branch (all, or those given as "from-to" or
        "from-to ckt") gets a fault at its from bus at fault_start,
        cleared after fault_duration by tripping it, replacing the base
        case's events, and runs dsf to end_time (s).  ranks is per job;
        the server runs --parallel jobs at a time.
        """
        entry = self._entry(case)
        if kind not in study.KINDS:
            raise ValueError("kind is one of %s, not %r"
                             % (", ".join(study.KINDS), kind))
        app = "powerflow" if kind == "load_sweep" else "dsf"
        if app not in entry["apps"]:
            raise ValueError("%s cannot run %s; it supports %s"
                             % (case, app, ", ".join(entry["apps"])))
        if not isinstance(ranks, int) or isinstance(ranks, bool) or ranks < 1:
            raise ValueError("ranks must be a positive int, not %r" % (ranks,))
        if ranks > entry["max_ranks"]:
            raise ValueError("%s (%d buses) runs on at most %d ranks"
                             % (case, entry["buses"], entry["max_ranks"]))
        if kind == "load_sweep":
            if not isinstance(factors, list) or not factors or any(
                    isinstance(f, bool) or not isinstance(f, (int, float))
                    or f <= 0 for f in factors):
                raise ValueError("factors must be a non-empty list of positive "
                                 "numbers, such as [1.0, 1.1, 1.2]")
            points = [{"name": "load-%g" % f, "factor": f,
                       "edits": [{"op": "scale_all_loads", "factor": f}] + (
                           [{"op": "scale_all_generation", "factor": f}]
                           if scale_generation else [])} for f in factors]
            params = {"factors": factors, "scale_generation": scale_generation}
            settings = {}
        else:
            c = self._load(case)
            raw = c.files()[next(p for p in c.files()
                                 if ".networkConfiguration" in p)]
            outages = study.n1_branches(raw, branches)
            end = end_time if end_time is not None else float(
                c.get("Dynamic_simulation.simulationTime"))
            if fault_start + fault_duration >= end:
                raise ValueError("the fault (%g s to %g s) must end before "
                                 "end_time %g s" % (fault_start, fault_start
                                                    + fault_duration, end))
            points = [{"name": "out-%d-%d-%s" % b, "branch": "%d-%d %s" % b,
                       "edits": [
                           {"op": "add_fault", "from": b[0], "to": b[1],
                            "ckt": b[2], "start": fault_start,
                            "duration": fault_duration},
                           {"op": "trip_branch", "from": b[0], "to": b[1],
                            "ckt": b[2], "time": fault_start + fault_duration}]}
                      for b in outages]
            params = {"fault_start": fault_start,
                      "fault_duration": fault_duration, "end_time": end,
                      "branches": [pt["branch"] for pt in points]}
            settings = {"end_time": end}
        width = min(self.runner.max_parallel, len(points))
        cpus = os.cpu_count() or 1
        if ranks * width > cpus:
            raise ValueError("%d jobs at a time on %d ranks each exceed the %d "
                             "cpus; lower ranks or start gridpack-mcp with a "
                             "smaller --parallel" % (width, ranks, cpus))
        st = study.Study(kind, case, params)
        # Every derived case is written before any job starts, so a bad
        # point fails the whole study with nothing running.
        for pt in points:
            d = self.modify_case(case, pt["edits"],
                                 name="%s-%s" % (st.id, pt["name"]),
                                 replace_events=(kind == "dynamic_n1"))
            pt["case"] = d["case"]
            if kind == "dynamic_n1":
                # Angle spread needs every machine, not the base's watch list.
                e = self._entry(d["case"])
                watch_all_generators(e["_xml"], raw)
        if kind == "dynamic_n1":
            with self._lock:
                self._cases = None            # re-read the rewritten XMLs
        for pt in points:
            job = self._submit(app, self._entry(pt["case"]), ranks, settings)
            st.items.append(dict({k: v for k, v in pt.items() if k != "edits"},
                                 job_id=job.id))
        with self._lock:
            self._studies[st.id] = st
        deadline = time.time() + max(0.0, min(float(wait_seconds), _MAX_WAIT))
        for it in st.items:
            self.runner.wait(it["job_id"], max(0.0, deadline - time.time()))
        return st.status(self.runner)

    def study_status(self, study_id: str) -> Dict[str, Any]:
        """Job states of a study, per item and in total."""
        return self._study(study_id).status(self.runner)

    def study_result(self, study_id: str) -> Dict[str, Any]:
        """The finished items ranked: load sweeps by factor with load,
        generation, losses, lowest voltage and violations; dynamic N-1 most
        severe first by loss of synchronism, peak rotor-angle spread and
        speed deviation.  Pending items are counted, failed ones listed."""
        return study.rank(self._study(study_id), self.runner)

    def _truth(self, case):
        """measurements.json of a synthesized SE case, or None."""
        if not case:
            return None
        try:
            path = Path(self._entry(case)["_xml"]).with_name("measurements.json")
        except ValueError:
            return None
        return json.loads(path.read_text()) if path.is_file() else None

    def _study(self, study_id):
        with self._lock:
            if study_id not in self._studies:
                raise KeyError("no study %r" % study_id)
            return self._studies[study_id]

    # ---- post-processing -------------------------------------------------

    def job_series(self, job_id: str, columns: Optional[List[str]] = None,
                   max_points: int = 200) -> Dict[str, Any]:
        """Chosen channels of a finished dynamic run (dsf watch, hadrec
        observations), decimated to at most max_points rows.

        columns: channel names or globs such as ``"*_speed"`` (time is
        always the first column); None means every channel.
        """
        job, result = self._finished(job_id)
        return post.series(result, job.app, columns,
                           max(2, min(int(max_points), 100000)))

    def job_report(self, job_id: str) -> str:
        """A short text verdict on a finished job: convergence, violations,
        totals, events, speed deviation, damping and angle spread."""
        job, result = self._finished(job_id)
        text = post.report(result, job.app, job.workdir)
        if job.app == "se":
            truth = self._truth(self._job_case.get(job_id))
            if truth:
                acc = measurements.accuracy(result, truth["truth_buses"],
                                            truth["planted"],
                                            result.get("bad_measurements", []))
                text = "\n".join([text] + measurements.estimate_lines(acc))
        return text

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
        if app == "hadrec" and not any(
                _records(self._load(case)._root)["observations"].values()):
            raise ValueError("%s has no observations, so hadrec would record "
                             "nothing; derive a case with modify_case "
                             "add_observations {\"generators\": \"all\", "
                             "\"buses\": \"all\"}" % case)
        job = self._submit(app, entry, ranks, settings)
        self.runner.wait(job.id, max(0.0, min(float(wait_seconds),
                                              _MAX_WAIT)))
        return self.job_status(job.id)

    def _submit(self, app, entry, ranks, settings):
        c = self._load(entry["case"])
        settings = {k: v for k, v in settings.items() if v is not None}
        if settings:
            c.override(app, **settings)
        # EMT aborts without its PETSc options file in the working dir.
        files = [entry["_petscrc"]] if app == "emt" and entry["_petscrc"] \
            else []
        job = self.runner.submit(app, c, ranks=ranks, files=files,
                                 timeout=self.timeout)
        with self._lock:
            self._job_case[job.id] = entry["case"]
        return job

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


def _records(root):
    """What a dynamic case records: watched generators and observations."""
    ds = root.find("Dynamic_simulation")
    watch = ds.find("generatorWatch") if ds is not None else None
    obs = ds.find("observations") if ds is not None else None
    counts = {}
    for o in (obs if obs is not None else []):
        kind = (o.findtext("type") or "").strip()
        counts[kind] = counts.get(kind, 0) + 1
    return {"watched_generators": len(watch) if watch is not None else 0,
            "observations": counts}


def _check_observe(e):
    extra = [k for k in e if k not in ("op", "generators", "buses", "bus_frequency")]
    if extra:
        raise ValueError("add_observations: unexpected %s" % ", ".join(extra))
    if e.get("generators") is None and e.get("buses") is None:
        raise ValueError("add_observations needs generators and/or buses "
                         "(\"all\" or a list of bus numbers)")
    for k in ("generators", "buses"):
        v = e.get(k)
        if v is not None and v != "all" and not (
                isinstance(v, list) and v and all(
                    isinstance(b, int) and not isinstance(b, bool) for b in v)):
            raise ValueError("add_observations: %s is \"all\" or a list of bus "
                             "numbers" % k)
    if not isinstance(e.get("bus_frequency", False), bool):
        raise ValueError("add_observations: bus_frequency is true or false")


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
        entry = {"case": name, "apps": apps,
                 "network": Path(net).name if net else None,
                 "buses": buses, "max_ranks": _max_ranks(buses),
                 "_xml": xml, "_base": base,
                 "_petscrc": local if local.is_file() else petscrc}
        if "emt" in apps:
            entry["notes"] = ["emt results are status only"]
        out.append(entry)
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
