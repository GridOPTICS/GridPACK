# -------------------------------------------------------------
# file: python/examples/_case.py
# -------------------------------------------------------------
# Shared helpers for the demos.
# -------------------------------------------------------------

"""Shared helpers for the demos."""

import contextlib
import csv
import ctypes
import os
import sys
import tempfile
from pathlib import Path


def stage(session, case, deferred=()) -> str:
    """Write `case` into this demo's directory, chdir there, return the XML name.

    Named after the case, not mkdtemp, so every rank gets the same path.
    """
    dest = Path(tempfile.gettempdir()) / ("gridpack-demo-" + case.source.stem)
    xml = case.write(dest, session=session, deferred=deferred)
    os.chdir(dest)
    return xml.name


def write_contingency_list(path: str, branches) -> None:
    """Write one single-branch outage per row of a branches() table.

    The format ca.x reads and ContingencyAnalysis parses; the circuit id
    is the branch's own, so multi-circuit corridors stay distinct.
    """
    lines = ['<?xml version="1.0" encoding="utf-8"?>', "<ContingencyList>",
             "  <Contingency_analysis>", "    <Contingencies>"]
    for br in branches:
        lines += [
            "      <Contingency>",
            "        <contingencyType>Line</contingencyType>",
            "        <contingencyName>%s</contingencyName>" % branch_name(br),
            "        <contingencyLineBuses> %d %d </contingencyLineBuses>"
            % (br["fromBus"], br["toBus"]),
            "        <contingencyLineNames> %s </contingencyLineNames>"
            % br["circuitId"].strip(),
            "      </Contingency>",
        ]
    lines += ["    </Contingencies>", "  </Contingency_analysis>",
              "</ContingencyList>", ""]
    Path(path).write_text("\n".join(lines))


def branch_name(br) -> str:
    """`L<from>-<to>_<ckt>`: the contingency name for a branch row."""
    return "L%d-%d_%s" % (br["fromBus"], br["toBus"], br["circuitId"].strip())


def write_csv(path: str, rows, columns=None) -> None:
    """Write already-gathered rows.  Call on one rank only.

    The results API's own to_csv() gathers *and* writes with no rank
    guard, so on N ranks it is N writers on one path; guarding the call
    instead deadlocks, because the gather is inside it.
    """
    rows = list(rows)
    columns = list(columns or (rows[0].keys() if rows else []))
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=columns)
        w.writeheader()
        w.writerows(rows)


@contextlib.contextmanager
def stdout_to(path: str):
    """Send everything written to fd 1 inside the block, C++ included, to
    `path`.  For modules that print with bare printf, which NoPrint and
    suppress_output cannot reach.  Both sides flush C stdio, or buffered
    output lands on the wrong side of the switch.
    """
    libc = ctypes.CDLL(None)
    sys.stdout.flush()
    libc.fflush(None)
    saved = os.dup(1)
    with open(path, "a") as log:
        os.dup2(log.fileno(), 1)
        try:
            yield
        finally:
            sys.stdout.flush()
            libc.fflush(None)
            os.dup2(saved, 1)
            os.close(saved)
