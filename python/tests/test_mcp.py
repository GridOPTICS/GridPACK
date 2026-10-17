# -------------------------------------------------------------
# file: python/tests/test_mcp.py
# -------------------------------------------------------------
# gridpack-mcp: install, version and the stdio handshake.
# -------------------------------------------------------------

import json
import shutil
import subprocess
import sys

import pytest

_SCRIPT = shutil.which("gridpack-mcp")
_CLI = pytest.mark.skipif(_SCRIPT is None,
                          reason="gridpack-mcp console script not installed")


@_CLI
def test_version_is_the_package_version():
    from importlib.metadata import version
    r = subprocess.run([_SCRIPT, "--version"], capture_output=True,
                       text=True, timeout=60)
    assert r.returncode == 0
    assert r.stdout.strip() == "gridpack-mcp %s" % version("gridpack")


@_CLI
def test_stdio_handshake():
    pytest.importorskip("mcp")
    p = subprocess.Popen([_SCRIPT], stdin=subprocess.PIPE,
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                         text=True)

    def send(msg):
        p.stdin.write(json.dumps(msg) + "\n")
        p.stdin.flush()

    try:
        send({"jsonrpc": "2.0", "id": 1, "method": "initialize",
              "params": {"protocolVersion": "2025-06-18", "capabilities": {},
                         "clientInfo": {"name": "pytest", "version": "0"}}})
        # Any stray print on stdout breaks this parse.
        init = json.loads(p.stdout.readline())
        send({"jsonrpc": "2.0", "method": "notifications/initialized"})
        send({"jsonrpc": "2.0", "id": 2, "method": "tools/list"})
        tools = json.loads(p.stdout.readline())
        p.stdin.close()
        assert p.wait(timeout=30) == 0
        for line in p.stdout.read().splitlines():
            json.loads(line)
    finally:
        p.kill()
    assert init["result"]["serverInfo"]["name"] == "gridpack"
    assert tools["result"]["tools"] == []


def test_missing_mcp_is_a_usage_error():
    r = subprocess.run(
        [sys.executable, "-c",
         "import sys; sys.modules['mcp'] = None\n"
         "from gridpack.mcp.server import main; sys.exit(main([]))"],
        capture_output=True, text=True, timeout=60)
    assert r.returncode == 2
    assert "pip install 'gridpack[mcp]'" in r.stderr
