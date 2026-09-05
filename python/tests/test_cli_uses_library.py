# -------------------------------------------------------------
# file: python/tests/test_cli_uses_library.py
# -------------------------------------------------------------
# Build guard: the CLI must drive gridpack.* wrappers, not the raw pybind11
# layer -- the library grew fixes the CLI's own copy never got.
# -------------------------------------------------------------

import ast
from pathlib import Path

import pytest

CLI_DIR = Path(__file__).resolve().parents[1] / "gridpack" / "cli"

# Matched on the tail of the dotted expression, so "hadrec.Module" catches
# gridpack.hadrec.Module.  CoarseTimer is deliberately absent: no wrapper,
# nothing to get wrong.
FORBIDDEN = (
    "Environment",
    "Communicator",
    "Configuration",
    "NoPrint",
    "TaskManager",
    "TaskCounter",
    "DSFullApp",
    "EventVector",
    "powerflow.Powerflow",
    "state_estimation.SEApp",
    "hadrec.Module",
    "emt.EMT",
)

# (function, dotted) pairs that may still appear.  cmd_dsf/cmd_hadrec are
# migrated in the next commit; only the cmd_emt entry should survive it.
ALLOWED = {
    ("cmd_dsf", "DSFullApp"),
    ("cmd_dsf", "gridpack.NoPrint"),
    ("cmd_hadrec", "gridpack.NoPrint"),
    ("cmd_hadrec", "gridpack.hadrec.Module"),
    ("cmd_hadrec", "gridpack.dynamic_simulation.EventVector"),
    ("cmd_emt", "gridpack.emt.EMT"),
}


def _is_forbidden(dotted: str) -> bool:
    if "_gridpack" in dotted:
        return True
    return any(dotted == f or dotted.endswith("." + f) for f in FORBIDDEN)


def _enclosing(tree, lineno: int) -> str:
    """Name of the innermost function containing ``lineno``, else "<module>"."""
    best = "<module>"
    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            if node.lineno <= lineno <= (node.end_lineno or node.lineno):
                best = node.name
    return best


def _low_level_uses(path: Path):
    """Every low-level name in ``path`` as (function, dotted, lineno).

    ast, not grep: main.py's own docstring names gridpack._gridpack to
    explain this rule and must not match.
    """
    tree = ast.parse(path.read_text())
    found = {}
    for node in ast.walk(tree):
        if isinstance(node, (ast.Name, ast.Attribute)):
            dotted = ast.unparse(node)
        elif isinstance(node, ast.ImportFrom):
            for alias in node.names:
                if _is_forbidden(alias.name):
                    found[(_enclosing(tree, node.lineno), alias.name)] = node.lineno
            continue
        elif isinstance(node, ast.Import):
            for alias in node.names:
                if _is_forbidden(alias.name):
                    found[(_enclosing(tree, node.lineno), alias.name)] = node.lineno
            continue
        else:
            continue
        if _is_forbidden(dotted):
            found[(_enclosing(tree, node.lineno), dotted)] = node.lineno
    return [(fn, dotted, line) for (fn, dotted), line in sorted(found.items())]


def _cli_files():
    files = sorted(CLI_DIR.glob("*.py"))
    assert files, f"no CLI sources under {CLI_DIR}"
    return files


@pytest.mark.parametrize("path", _cli_files(), ids=lambda p: p.name)
def test_cli_does_not_use_low_level_api(path):
    """Fail the build when the CLI reaches past the wrapper classes."""
    bad = [(fn, dotted, line) for fn, dotted, line in _low_level_uses(path)
           if (fn, dotted) not in ALLOWED]
    assert not bad, "\n".join(
        "%s:%d: %s uses %s -- add a gridpack.* wrapper instead"
        % (path.name, line, fn, dotted) for fn, dotted, line in bad)


def test_allowlist_has_no_stale_entries():
    """A migrated command must shrink ALLOWED, or the guard rots open."""
    live = {(fn, dotted) for path in _cli_files()
            for fn, dotted, _ in _low_level_uses(path)}
    assert not (ALLOWED - live), (
        "ALLOWED lists names the CLI no longer uses: %s" % sorted(ALLOWED - live))
