# -------------------------------------------------------------
# file: gridpack/mcp/server.py
# -------------------------------------------------------------
# gridpack-mcp entry point (stdio transport).
# -------------------------------------------------------------

import argparse
import functools
import math
import sys

from ..cli.main import _package_version

_INSTALL_HINT = ("gridpack-mcp needs the mcp package: "
                 "pip install 'gridpack[mcp]' (Python 3.10+)\n")


_TOOLS = ("list_cases", "describe_case", "create_case", "create_dyr",
          "modify_case", "create_measurements", "run_powerflow", "run_contingency", "run_dynamic_sim",
          "run_state_estimation", "run_hadrec", "run_emt", "job_status",
          "job_result", "job_log", "job_series", "job_report", "compare_jobs",
          "plot_job", "run_study", "study_status", "study_result",
          "save_figure", "save_report", "export_data")


def build_server(tools=None):
    """The MCP server object, with ``tools`` (a :class:`Tools`) and the
    study prompts registered."""
    from mcp.server.mcpserver import MCPServer
    from .tools import Tools
    tools = tools or Tools()
    server = MCPServer(
        name="gridpack", version=_package_version(),
        instructions="Run GridPACK on its cases: list_cases (with each "
                     "case's max_ranks), then a run_* tool, then "
                     "job_status / job_result.  create_case makes a case "
                     "from a RAW (and DYR) file under a --cases directory, "
                     "create_dyr a generic DYR for a RAW; "
                     "modify_case derives one with loads, generators, "
                     "shunts or branches changed and faults or trips added.  "
                     "job_report, job_series, plot_job and compare_jobs read "
                     "finished jobs back; run_study fans out a load sweep or "
                     "a dynamic N-1 screen; create_measurements turns a power "
                     "flow into a state estimation case; save_figure, "
                     "save_report and export_data write files for papers and "
                     "readers.  EMT results are status only.")
    for name in _TOOLS:
        fn = getattr(tools, name)
        server.add_tool(_reported(_as_image(fn) if name == "plot_job" else fn))
    _prompts(server)
    return server


_WORKFLOW = """Work with the gridpack tools in this order and report what each
step returned.
1. list_cases, then describe_case for the case: note max_ranks (never exceed
   it; a 9-bus case runs on at most 3 ranks) and, for dynamics, records
   (watched generators, observations) and settings (time step, end time).
2. If the user has a RAW/DYR file under a --cases directory, create_case first.
3. Change the network or add events with modify_case; use replace_events
   when the base case already has a fault.  Event times must be multiples of
   the time step.  HADREC needs a fault as the first event and observations
   (add_observations).
4. Run with run_powerflow / run_dynamic_sim / run_hadrec / run_contingency;
   poll job_status until done when a run waits longer than 10 s.
5. Read results with job_report first, then plot_job or job_series, and
   compare_jobs against the base run.  Derived cases disappear when the
   server restarts; each edits.json beside a derived case says how it was made.
"""


def _prompts(server):
    @server.prompt(name="stability_study",
                   description="Fault a branch of a case, run the dynamic "
                               "simulation and judge stability")
    def stability_study(case: str, branch: str = "",
                        fault_start: str = "0.1",
                        fault_duration: str = "0.05") -> str:
        where = ("branch %s" % branch) if branch else "a branch you choose"
        return (_WORKFLOW + "\nTask: on case %s, put a fault on %s at %s s "
                "for %s s cleared by tripping the branch (modify_case with "
                "add_fault then trip_branch, replace_events true), run the "
                "dynamic simulation, give job_report's verdict, plot the "
                "generator speeds, and compare with the unfaulted run."
                % (case, where, fault_start, fault_duration))

    @server.prompt(name="n1_screen",
                   description="Screen every branch outage of a case with "
                               "dynamic simulations and rank them")
    def n1_screen(case: str, end_time: str = "") -> str:
        return (_WORKFLOW + "\nTask: run_study kind dynamic_n1 on %s%s, poll "
                "study_status, then study_result; list the outages ranked "
                "most severe first, explain unstable (angle spread beyond "
                "180 deg or a machine still running away at the end) as a "
                "first-swing screen, and plot the worst one."
                % (case, (" to %s s" % end_time) if end_time else ""))

    @server.prompt(name="voltage_study",
                   description="Sweep load growth on a case and report "
                               "voltages, losses and violations")
    def voltage_study(case: str, growth_percent: str = "30") -> str:
        return (_WORKFLOW + "\nTask: run_study kind load_sweep on %s with "
                "factors from 1.0 up to 1 + %s/100 in steps of 0.05 (or "
                "coarser), then study_result; report the lowest voltage and "
                "its bus, violations and losses per factor, and the factor "
                "at which the first violation or divergence appears."
                % (case, growth_percent))


def _as_image(plot_job):
    """plot_job returning the PNG as MCP image content."""
    from typing import List, Optional
    from mcp.server.mcpserver.utilities.types import Image

    def tool(job_id: str, kind: Optional[str] = None,
             channels: Optional[List[str]] = None,
             title: Optional[str] = None) -> Image:
        return Image(data=plot_job(job_id, kind, channels, title), format="png")
    tool.__name__, tool.__doc__ = "plot_job", plot_job.__doc__
    return tool


def _reported(fn):
    """Pass expected failures to the client; mcp masks other exceptions."""
    from mcp.server.mcpserver.exceptions import ToolError

    @functools.wraps(fn)
    def tool(*args, **kwargs):
        try:
            return fn(*args, **kwargs)
        except (ValueError, KeyError, TypeError) as exc:
            raise ToolError(str(exc.args[0] if exc.args else exc)) from exc
    return tool


def _seconds(text):
    value = float(text)
    if not math.isfinite(value) or value <= 0:
        raise argparse.ArgumentTypeError("must be a positive number of seconds")
    return value


def _count(text):
    value = int(text)
    if value < 1:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return value


def _doctor(args):
    """Check the installation; return 0 when the server can serve."""
    import os
    import platform
    import shutil
    import subprocess
    ok = True

    def line(good, text):
        print(("ok       " if good else "MISSING  ") + text)
        return good

    ok &= line(sys.version_info >= (3, 10), "Python %s (3.10+ needed)"
               % platform.python_version())
    try:
        import mcp  # noqa: F401
        line(True, "mcp package")
    except ImportError:
        ok = line(False, "mcp package: pip install 'gridpack[mcp]'")
    cli = shutil.which("gridpack")
    if cli:
        v = subprocess.run([cli, "--version"], capture_output=True, text=True)
        ok &= line(v.returncode == 0, "gridpack CLI %s (%s)"
                   % (cli, (v.stdout or v.stderr).strip()))
    else:
        ok = line(False, "gridpack CLI on PATH: the venv's bin directory must "
                         "be on the server's PATH")
    mpi = shutil.which("mpiexec")
    line(bool(mpi), "mpiexec %s" % (mpi or "on PATH: runs stay on one rank"))
    try:
        import matplotlib  # noqa: F401
        line(True, "matplotlib (plot_job)")
    except ImportError:
        line(False, "matplotlib: pip install 'gridpack[plot]' (plot_job refuses "
                    "without it)")
    from .runner import JobRunner
    from .tools import Tools
    try:
        tools = Tools(JobRunner(max_parallel=args.parallel), cases=args.cases,
                      output=args.output)
        cases = tools.list_cases()
        own = {}
        for c in cases:
            label = c["case"].split("/")[0]
            if label not in ("input", "emt"):
                own[label] = own.get(label, 0) + 1
        line(True, "%d bundled cases%s" % (
            len(cases) - sum(own.values()),
            "".join("; %d under --cases %s" % (n, d) for d, n in own.items())))
    except ValueError as exc:
        ok = line(False, str(exc))
    cpus = os.cpu_count() or 1
    line(True, "%d cpus: --parallel %d jobs x ranks per job must stay within"
         % (cpus, args.parallel))
    line(True, "output: %s" % (args.output or "temporary (use --output DIR to keep "
                                              "figures and reports)"))
    print("ready" if ok else "not ready")
    return 0 if ok else 1


def _config(args):
    """Print how to register this server with Claude Code and Claude Desktop."""
    import json
    import os
    import shutil
    cmd = shutil.which("gridpack-mcp") or os.path.abspath(sys.argv[0])
    opts = []
    if args.parallel != 1:
        opts += ["--parallel", str(args.parallel)]
    if args.timeout != 600.0:
        opts += ["--timeout", "%g" % args.timeout]
    for d in args.cases:
        opts += ["--cases", os.path.abspath(os.path.expanduser(d))]
    if args.output:
        opts += ["--output", os.path.abspath(os.path.expanduser(args.output))]
    dirs = []
    for exe in (cmd, shutil.which("gridpack"), shutil.which("mpiexec")):
        d = os.path.dirname(exe) if exe else None
        if d and d not in dirs:
            dirs.append(d)
    path = ":".join(dirs + ["/usr/bin", "/bin"])
    print("# Claude Code")
    print("claude mcp add gridpack -e PATH=%s -- %s" % (path, " ".join([cmd] + opts)))
    print("\n# Claude Desktop: merge into claude_desktop_config.json")
    print(json.dumps({"mcpServers": {"gridpack": {
        "command": cmd, "args": opts, "env": {"PATH": path}}}}, indent=2))
    return 0


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="gridpack-mcp",
        description="GridPACK MCP server over stdio.  With no command it "
                    "serves; 'doctor' checks the installation; 'config' "
                    "prints the client registration for the given options.")
    parser.add_argument("command", nargs="?", choices=("doctor", "config"),
                        help="doctor | config (default: serve)")
    parser.add_argument("--version", action="version",
                        version="gridpack-mcp %s" % _package_version())
    parser.add_argument("--timeout", type=_seconds, default=600.0,
                        metavar="SECONDS",
                        help="kill a job after this long (default 600)")
    parser.add_argument("--parallel", type=_count, default=1, metavar="N",
                        help="run up to N jobs at a time (default 1)")
    parser.add_argument("--output", metavar="DIR",
                        help="where save_figure, save_report and export_data "
                             "write (default: the server's temporary directory)")
    parser.add_argument("--cases", action="append", default=[],
                        metavar="DIR",
                        help="also offer the XMLs under DIR as cases, named "
                             "<DIR name>/<path>; repeatable")
    args = parser.parse_args(argv)
    if args.command == "doctor":
        return _doctor(args)
    if args.command == "config":
        return _config(args)
    from .runner import JobRunner
    from .tools import Tools
    try:
        tools = Tools(JobRunner(max_parallel=args.parallel),
                      timeout=args.timeout, cases=args.cases, output=args.output)
        tools.list_cases()          # report skipped cases at startup
    except ValueError as exc:
        sys.stderr.write("gridpack-mcp: %s\n" % exc)
        return 2
    try:
        server = build_server(tools)
    except ImportError:
        sys.stderr.write(_INSTALL_HINT)
        return 2
    try:
        # stdout carries the protocol: nothing else may print to it.
        server.run("stdio")
    finally:
        # Jobs run in their own sessions and would outlive the server.
        for job in tools.runner.jobs():
            tools.runner.cancel(job.id)
    return 0


if __name__ == "__main__":
    sys.exit(main())
