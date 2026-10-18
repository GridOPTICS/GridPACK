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


_TOOLS = ("list_cases", "describe_case", "run_powerflow", "run_contingency",
          "run_dynamic_sim", "run_state_estimation", "run_hadrec", "run_emt",
          "job_status", "job_result", "job_log")


def build_server(tools=None):
    """The MCP server object, with ``tools`` (a :class:`Tools`) registered."""
    from mcp.server.mcpserver import MCPServer
    from .tools import Tools
    tools = tools or Tools()
    server = MCPServer(
        name="gridpack", version=_package_version(),
        instructions="Run GridPACK on its bundled cases: list_cases (with "
                     "each case's max_ranks), then a run_* tool, then "
                     "job_status / job_result.")
    for name in _TOOLS:
        server.add_tool(_reported(getattr(tools, name)))
    return server


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


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="gridpack-mcp",
        description="GridPACK MCP server over stdio.")
    parser.add_argument("--version", action="version",
                        version="gridpack-mcp %s" % _package_version())
    parser.add_argument("--timeout", type=_seconds, default=600.0,
                        metavar="SECONDS",
                        help="kill a job after this long (default 600)")
    args = parser.parse_args(argv)
    from .tools import Tools
    tools = Tools(timeout=args.timeout)
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
