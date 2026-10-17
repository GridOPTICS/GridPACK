# -------------------------------------------------------------
# file: gridpack/mcp/server.py
# -------------------------------------------------------------
# gridpack-mcp entry point (stdio transport).
# -------------------------------------------------------------

import argparse
import sys

from ..cli.main import _package_version

_INSTALL_HINT = ("gridpack-mcp needs the mcp package: "
                 "pip install 'gridpack[mcp]' (Python 3.10+)\n")


def build_server():
    """The MCP server object."""
    from mcp.server.mcpserver import MCPServer
    return MCPServer(name="gridpack", version=_package_version())


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
