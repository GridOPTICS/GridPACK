# -------------------------------------------------------------
# file: gridpack/mcp/server.py
# -------------------------------------------------------------
# gridpack-mcp entry point (stdio transport).
# -------------------------------------------------------------

import argparse
import functools
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
        instructions="Run GridPACK on its bundled cases: list_cases, then a "
                     "run_* tool, then job_status / job_result.")
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


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="gridpack-mcp",
        description="GridPACK MCP server over stdio.")
    parser.add_argument("--version", action="version",
                        version="gridpack-mcp %s" % _package_version())
    parser.parse_args(argv)
    try:
        server = build_server()
    except ImportError:
        sys.stderr.write(_INSTALL_HINT)
        return 2
    # stdout carries the protocol: nothing else may print to it.
    server.run("stdio")
    return 0


if __name__ == "__main__":
    sys.exit(main())
