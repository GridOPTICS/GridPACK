# -------------------------------------------------------------
# file: python/examples/_case.py
# -------------------------------------------------------------
# Shared input staging for the demos.
# -------------------------------------------------------------

"""Copy a demo's inputs into one directory and work there.

GridPACK resolves ``<networkConfiguration>`` relative to the working
directory, and the stock inputs keep the XML and the network file in
separate trees, so the files have to be brought together first.  A
missing network file aborts with ``wrong dimension specified``, not a
file-not-found error, which is hard to place if you have not seen it.
"""

import csv
import os
import shutil
import tempfile
from pathlib import Path


def data_sets() -> Path:
    """The GridPACK source tree's data_sets directory.

    ``GRIDPACK_TEST_DATA`` overrides the search for installs without the
    source tree, matching what the test suite honours.
    """
    env = os.environ.get("GRIDPACK_TEST_DATA")
    root = Path(env).resolve() if env else Path(__file__).resolve().parents[2]
    d = root / "src/applications/data_sets"
    if not d.is_dir():
        raise SystemExit(
            "no input data at %s\n"
            "Set GRIDPACK_TEST_DATA to the GridPACK source tree." % d)
    return d


def stage(session, *relpaths: str) -> str:
    """Stage inputs into a shared directory, chdir there, return the XML name.

    One writer and a barrier: every rank chdirs to the same path, so the
    directory name is derived from the case rather than from mkdtemp,
    which would hand each rank a different one.
    """
    dest = Path(tempfile.gettempdir()) / ("gridpack-demo-" + Path(relpaths[0]).stem)

    if session.rank == 0:
        src_dir = data_sets()
        dest.mkdir(parents=True, exist_ok=True)
        for rel in relpaths:
            src = src_dir / rel
            if not src.is_file():
                raise SystemExit("missing input: %s" % src)
            shutil.copy(src, dest / Path(rel).name)
    session.barrier()

    os.chdir(dest)
    return Path(relpaths[0]).name


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
