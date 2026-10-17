# -------------------------------------------------------------
# file: gridpack/case.py
# -------------------------------------------------------------
# Build, edit and stage an input XML from Python.
# -------------------------------------------------------------

"""Assemble one input file per study and stage it with its data files.

GridPACK keeps a single process-wide ``Configuration`` tree: every
``open()`` merges another root into it and a lookup searches only the
first, so whatever a process needs must be in the first file it opens.
A :class:`Case` holds that one file in memory until :meth:`Case.write`.
"""

import copy
import os
import re
import shutil
import xml.etree.ElementTree as ET
from pathlib import Path

__all__ = ["Case", "data_sets"]

# Tags naming an input file, and the data_sets directory that holds it.
_INPUT_TAGS = (
    (re.compile(r"networkConfiguration(_v\d+)?$"), "raw"),
    (re.compile(r"generatorParameters$"), "dyr"),
    (re.compile(r"contingencyList$"), "contingencies"),
    (re.compile(r"measurementList$"), "measurements"),
    (re.compile(r"Kalman(Ang|Mag)Data$"), "kalman"),
)


def data_sets() -> Path:
    """The ``src/applications/data_sets`` directory of a GridPACK tree.

    Searched in ``GRIDPACK_TEST_DATA``, then ``GRIDPACK_ROOT_DIR`` (set in
    the Docker image), then the source tree this package was installed from.
    """
    roots = [os.environ.get("GRIDPACK_TEST_DATA"),
             os.environ.get("GRIDPACK_ROOT_DIR"),
             Path(__file__).resolve().parents[2]]
    for root in roots:
        if root and (Path(root) / "src/applications/data_sets").is_dir():
            return (Path(root) / "src/applications/data_sets").resolve()
    raise FileNotFoundError(
        "no GridPACK data_sets directory found; set GRIDPACK_TEST_DATA to "
        "the GridPACK source tree")


def _input_dir(tag):
    for pattern, subdir in _INPUT_TAGS:
        if pattern.match(tag):
            return subdir
    return None


def _parse(text, origin):
    """Parse what Boost's XML reader accepts.

    It skips the declaration and any ``<!...>`` node, so the stock
    ``<!-ksp_view>`` and a misspelt declaration are harmless to GridPACK
    and fatal to ElementTree.  Dropping them gives the tree GridPACK sees.
    """
    text = re.sub(r"<\?.*?\?>", "", text, flags=re.S)
    text = re.sub(r"<!--.*?-->", "", text, flags=re.S)
    text = re.sub(r"<!(?!\[CDATA\[)[^>]*>", "", text)
    try:
        root = ET.fromstring(text)
    except ET.ParseError as exc:
        raise ValueError("%s: %s" % (origin, exc)) from None
    if root.tag != "Configuration":
        raise ValueError("%s: root is <%s>, not <Configuration>"
                         % (origin, root.tag))
    return root


def _format(value):
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float):
        return repr(value)
    return str(value)


class Case:
    """An input XML held in memory, edited by dotted path, written once.

    Usage::

        case = gridpack.Case.from_data_sets("input/powerflow/input_14.xml")
        case.set("Powerflow.tolerance", 1e-8)
        xml = case.write("/tmp/run")      # input_14.xml + IEEE14.raw
        os.chdir(xml.parent)              # GridPACK reads relative to cwd

    Paths are relative to ``<Configuration>`` and name the first match at
    each level, as GridPACK's own lookup does.  Comments are not kept.

    File references (``networkConfiguration``, ``generatorParameters``,
    ``contingencyList``, ``measurementList``, Kalman data) are found
    next to the XML, then in the working directory, then in the matching
    ``data_sets`` subdirectory for a data_sets case.  :meth:`get` returns
    the resolved path; :meth:`write` copies each file beside the XML.
    """

    def __init__(self, path, *, _data_sets=None):
        path = Path(path).resolve()
        self.source = path
        self._root = _parse(path.read_text(encoding="utf-8"), path)
        self._base = path.parent
        self._data_sets = _data_sets
        for el in self._root.iter():
            if _input_dir(el.tag) and len(el) == 0:
                found = self._resolve(el.tag, el.text)
                if found:
                    el.text = str(found)

    @classmethod
    def from_data_sets(cls, relpath, root=None):
        """Load ``relpath`` from :func:`data_sets`, or from ``root``."""
        d = Path(root) if root else data_sets()
        return cls(d / relpath, _data_sets=d)

    # ---- edits -------------------------------------------------------

    def get(self, path, default=None):
        """Text at ``path``, stripped, or ``default`` when absent."""
        el = self._find(path)
        if el is None or el.text is None:
            return default
        return el.text.strip()

    def set(self, path, value):
        """Set the text at ``path``, creating missing elements.

        Booleans are written ``true``/``false``.  A file reference is
        resolved now; one not found yet must be named in
        :meth:`write`'s ``deferred``.
        """
        el = self._find(path, create=True)
        if len(el):
            raise ValueError("<%s> holds elements; use set_xml" % el.tag)
        text = _format(value)
        if _input_dir(el.tag):
            found = self._resolve(el.tag, text)
            text = str(found) if found else text
        el.text = text

    def set_xml(self, path, fragment):
        """Replace the content of ``path`` with an XML ``fragment``."""
        el = self._find(path, create=True)
        new = _parse("<Configuration>%s</Configuration>" % fragment,
                     "fragment for %s" % path)
        el.clear()
        el.extend(new)

    def include(self, other, block):
        """Copy the top-level ``block`` of ``other``, replacing any here.

        For composing one study from several stock inputs, whose blocks
        must all be in the one file opened first.  ``other`` is a Case or
        a data_sets path; its file references stay resolved.
        """
        if not isinstance(other, Case):
            other = Case.from_data_sets(other, root=self._data_sets)
        src = other._root.find(block)
        if src is None:
            raise KeyError("no <%s> in %s" % (block, other.source))
        old = self._root.find(block)
        if old is not None:
            self._root.remove(old)
        self._root.append(copy.deepcopy(src))

    def files(self):
        """``{path: value}`` for every file reference."""
        return {p: (el.text or "").strip() for p, el in self._references()}

    # ---- output ------------------------------------------------------

    def write(self, workdir, name=None, *, session=None, deferred=()):
        """Write the XML and its input files into ``workdir``.

        Returns the XML's path.  References are rewritten to bare names,
        so GridPACK must run with ``workdir`` as its working directory.
        ``deferred`` names files the caller puts in ``workdir`` before
        GridPACK reads them.  With ``session``, rank 0 writes and all
        ranks wait for it.
        """
        workdir = Path(workdir).resolve()
        root = copy.deepcopy(self._root)
        copies = {}
        for path, el in Case._references_in(root):
            value = (el.text or "").strip()
            src = self._resolve(el.tag, value)
            if src is None:
                if Path(value).name in deferred:
                    el.text = Path(value).name
                    continue
                raise FileNotFoundError(
                    "%s names %r, which is not next to %s, in the working "
                    "directory%s" % (path, value, self.source,
                                     " or in data_sets" if self._data_sets
                                     else ""))
            if copies.get(src.name, src) != src:
                raise ValueError("two different files named %s: %s and %s"
                                 % (src.name, copies[src.name], src))
            copies[src.name] = src
            el.text = src.name

        xml = workdir / (name or self.source.name)
        if xml.name in copies:
            raise ValueError("%s would overwrite an input file" % xml.name)
        if session is None or session.rank == 0:
            workdir.mkdir(parents=True, exist_ok=True)
            for dest, src in copies.items():
                if src != workdir / dest:
                    shutil.copyfile(src, workdir / dest)
            ET.indent(root, "  ")
            ET.ElementTree(root).write(xml, encoding="utf-8",
                                       xml_declaration=True)
        if session is not None:
            session.barrier()
        return xml

    # ---- internals -------------------------------------------------------

    def _find(self, path, create=False):
        parts = path.split(".")
        if parts[0] == "Configuration":
            parts = parts[1:]
        if not parts or not all(parts):
            raise ValueError("bad path %r" % path)
        el = self._root
        for tag in parts:
            child = el.find(tag)
            if child is None:
                if not create:
                    return None
                child = ET.SubElement(el, tag)
            el = child
        return el

    def _references(self):
        return Case._references_in(self._root)

    @staticmethod
    def _references_in(root):
        def walk(el, prefix):
            for child in el:
                path = prefix + child.tag
                if len(child):
                    yield from walk(child, path + ".")
                elif _input_dir(child.tag):
                    yield path, child
        return list(walk(root, ""))

    def _resolve(self, tag, value):
        value = (value or "").strip()
        if not value:
            return None
        p = Path(value)
        candidates = [p] if p.is_absolute() else [self._base / p, Path.cwd() / p]
        if self._data_sets and not p.is_absolute():
            candidates.append(self._data_sets / _input_dir(tag) / p.name)
        for c in candidates:
            if c.is_file():
                return c.resolve()
        return None
