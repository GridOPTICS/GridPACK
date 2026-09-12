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
import re
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


def stage(session, *relpaths: str, network: str = None) -> str:
    """Stage inputs into a shared directory, chdir there, return the XML name.

    One writer and a barrier: every rank chdirs to the same path, so the
    directory name is derived from the case rather than from mkdtemp,
    which would hand each rank a different one.

    ``network`` stages one more file and repoints the XML at it, for
    running a stock config against a different network than it names.
    """
    dest = Path(tempfile.gettempdir()) / ("gridpack-demo-" + Path(relpaths[0]).stem)
    wanted = list(relpaths) + ([network] if network else [])

    if session.rank == 0:
        src_dir = data_sets()
        dest.mkdir(parents=True, exist_ok=True)
        for rel in wanted:
            src = src_dir / rel
            if not src.is_file():
                raise SystemExit("missing input: %s" % src)
            shutil.copy(src, dest / Path(rel).name)
        if network:
            _repoint(dest / Path(relpaths[0]).name, Path(network).name)
    session.barrier()

    os.chdir(dest)
    return Path(relpaths[0]).name


def _repoint(xml: Path, network: str) -> None:
    """Rewrite the XML's <networkConfiguration> to name `network`."""
    set_tag(xml, "networkConfiguration", network)


def copy_input(session, relpath: str, name: str) -> str:
    """Stage one more data_sets file into the cwd under `name`.

    For a second XML whose stock filename collides with the first one's.
    One writer and a barrier, like stage().
    """
    if session.rank == 0:
        src = data_sets() / relpath
        if not src.is_file():
            raise SystemExit("missing input: %s" % src)
        shutil.copy(src, name)
    session.barrier()
    return name


def set_tag(xml, tag: str, value, parent: str = None) -> None:
    """Set the text of the first <tag> in `xml`, inserting it if absent.

    Comments are dropped first: the stock inputs keep alternatives such as
    <contingencyList> inside <!-- -->, and a plain substitution would edit
    the comment and leave the live config unchanged.  `parent` says where
    to insert when the tag is missing.
    """
    xml = Path(xml)
    text = re.sub(r"<!--.*?-->", "", xml.read_text(), flags=re.S)
    text, n = re.subn(r"(?<=<%s>).*?(?=</%s>)" % (tag, tag),
                      " %s " % value, text, count=1, flags=re.S)
    if not n:
        if parent is None:
            raise SystemExit("no <%s> in %s" % (tag, xml))
        text, n = re.subn(r"(<%s>)" % parent,
                          r"\1\n    <%s>%s</%s>" % (tag, value, tag),
                          text, count=1)
        if not n:
            raise SystemExit("no <%s> in %s" % (parent, xml))
    xml.write_text(text)


def copy_block(src_relpath: str, xml, tag: str) -> None:
    """Append <tag>...</tag> from a data_sets XML to the end of `xml`.

    For composing one config from two stock ones.  GridPACK holds a single
    process-wide Configuration tree and a path lookup only searches the
    first <Configuration> root, so a block an application needs has to be
    in the first file opened; a second file cannot supply it.
    """
    src = (data_sets() / src_relpath).read_text()
    m = re.search(r"<%s>.*?</%s>" % (tag, tag), src, flags=re.S)
    if not m:
        raise SystemExit("no <%s> in %s" % (tag, src_relpath))
    xml = Path(xml)
    text, n = re.subn(r"</Configuration>", "  %s\n</Configuration>" % m.group(0),
                      xml.read_text(), count=1)
    if not n:
        raise SystemExit("no </Configuration> in %s" % xml)
    xml.write_text(text)


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
