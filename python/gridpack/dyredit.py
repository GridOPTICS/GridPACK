# -------------------------------------------------------------
# file: gridpack/dyredit.py
# -------------------------------------------------------------
# Read, synthesize and edit PSS/E DYR dynamic model records.
# -------------------------------------------------------------

"""DYR records: ``bus, 'MODEL', 'id', p1, p2, ... /``.

A record runs to its ``/``; GridPACK splits it on commas and blanks.  The
field names below follow GridPACK's parsers (``src/parser/parser_classes``),
so an edit by name lands where the solver reads it.  Synthesized models use
generic textbook parameters: good for teaching, screening and comparing
scenarios, not validated for any real machine.
"""

import re
from pathlib import Path

from .rawedit import EditError, _id

__all__ = ["FIELDS", "MODEL_SETS", "dyr_records", "edit_dyr", "synthesize_dyr"]

# Parameter names in record order after bus, model and id (token 3 on).
FIELDS = {
    "GENCLS": ["H", "D", "XDP"],
    "GENROU": ["TDOP", "TDOPP", "TQOP", "TQOPP", "H", "D", "XD", "XQ", "XDP",
               "XQP", "XDPP", "XL", "S1", "S12"],
    "SEXS": ["TA_TB", "TB", "K", "TE", "EMIN", "EMAX"],
    "TGOV1": ["R", "T1", "VMAX", "VMIN", "T2", "T3", "DT"],
}

# Generic parameters on the machine's own MVA base (H in seconds).
MODEL_SETS = {
    "classical": [("GENCLS", [5.0, 2.0, 0.3])],
    "detailed": [("GENROU", [6.0, 0.05, 1.0, 0.05, 4.0, 0.0, 1.8, 1.7, 0.3,
                             0.55, 0.25, 0.15, 0.05, 0.3]),
                 ("SEXS", [0.1, 10.0, 100.0, 0.1, 0.0, 4.0]),
                 ("TGOV1", [0.05, 0.5, 1.05, 0.0, 2.1, 7.0, 0.0])],
}


def _tokens(body):
    """Fields split on commas and blanks outside quotes, quotes kept."""
    out, cur, quoted = [], [], False
    for ch in body:
        if ch == "'":
            quoted = not quoted
            cur.append(ch)
        elif (ch == "," or ch.isspace()) and not quoted:
            if cur:
                out.append("".join(cur))
                cur = []
        else:
            cur.append(ch)
    if cur:
        out.append("".join(cur))
    return out


class Record:
    """One DYR record: its tokens and where it sits in the file."""

    def __init__(self, text, start, end):
        self.text, self.start, self.end = text, start, end
        raw = _tokens(text.split("/", 1)[0])
        self.tokens = raw
        self.bus = int(raw[0])
        self.model = raw[1].strip("'\"").strip().upper()
        self.id = _id(raw[2])

    def names(self):
        """Field names: the parser's for known models, then ``p<token>``
        by position (p3 is the first parameter)."""
        known = FIELDS.get(self.model, [])
        return known + ["p%d" % i for i in range(len(known) + 3, len(self.tokens))]

    def get(self, field):
        return float(self.tokens[self._index(field)])

    def set(self, field, value):
        i = self._index(field)
        old = self.tokens[i]
        text = ("%.6f" % value).rstrip("0")
        text += "0" if text.endswith(".") else ""
        # Replace the token in the record text, keeping everything else.
        head = self.text.split("/", 1)
        pattern = re.compile(r"(?<![\w.\-])" + re.escape(old) + r"(?![\w.])")
        hits = list(pattern.finditer(head[0]))
        n = sum(1 for t in self.tokens[:i] if t == old)
        if len(hits) <= n:
            raise EditError("cannot place field %s in record %r" % (field, self.text))
        m = hits[n]
        head[0] = head[0][:m.start()] + text.rjust(len(old)) + head[0][m.end():]
        self.text = "/".join(head)
        self.tokens[i] = text

    def _index(self, field):
        names = self.names()
        if field not in names:
            raise EditError("%s has no field %s; fields are %s" % (
                self.model, field, ", ".join(names) or "none"))
        i = names.index(field) + 3
        if i >= len(self.tokens):
            raise EditError("record for %s %d id %s has no %s value" % (
                self.model, self.bus, self.id, field))
        return i


def dyr_records(path):
    """Every record in a DYR file, in order."""
    text = Path(path).read_text(errors="replace")
    out, pos = [], 0
    while True:
        end = text.find("/", pos)
        if end < 0:
            break
        chunk = text[pos:end + 1]
        # Keep trailing comment up to the line end with the record.
        nl = text.find("\n", end)
        tail_end = len(text) if nl < 0 else nl
        chunk = text[pos:tail_end]
        if chunk.strip() and chunk.split("/", 1)[0].strip():
            try:
                out.append(Record(chunk, pos, tail_end))
            except (IndexError, ValueError):
                raise EditError("unreadable DYR record: %r" % chunk.strip()[:60])
        pos = tail_end
    return out


def synthesize_dyr(generators, model_set="classical"):
    """DYR text giving every ``(bus, id)`` the models of ``model_set``."""
    if model_set not in MODEL_SETS:
        raise EditError("model_set is one of %s, not %r"
                        % (", ".join(MODEL_SETS), model_set))
    # No comment line: GridPACK's parser aborts on a record without data.
    lines = []
    for bus, gid in generators:
        for model, params in MODEL_SETS[model_set]:
            lines.append("%d, '%s', '%-2s', %s /" % (
                bus, model, gid, ", ".join("%g" % p for p in params)))
    return "\n".join(lines) + "\n"


def edit_dyr(src, dst, edits):
    """Apply ``[{bus, id, model, field, value}]`` to the DYR ``src`` and
    write ``dst``; return each edit with ``was`` and ``now``.  Every edit
    is checked before anything is written."""
    if not isinstance(edits, list) or not edits:
        raise EditError("edits must be a non-empty list")
    records = dyr_records(src)
    plan = []
    for n, e in enumerate(edits, 1):
        for k in ("bus", "model", "field", "value"):
            if k not in e:
                raise EditError("dyr edit %d: missing %s" % (n, k))
        if isinstance(e["value"], bool) or not isinstance(e["value"], (int, float)):
            raise EditError("dyr edit %d: value must be a number" % n)
        gid = _id(str(e.get("id", "1")))
        model = str(e["model"]).upper()
        hits = [r for r in records if r.bus == e["bus"] and r.id == gid
                and r.model == model]
        if not hits:
            have = sorted({r.model for r in records if r.bus == e["bus"] and r.id == gid})
            raise EditError("dyr edit %d: no %s record for generator %d id %s%s" % (
                n, model, e["bus"], gid,
                "; it has %s" % ", ".join(have) if have else ""))
        rec = hits[0]
        rec._index(e["field"])
        plan.append((e, rec))
    text = Path(src).read_text(errors="replace")
    report = []
    for e, rec in plan:
        was = rec.get(e["field"])
        rec.set(e["field"], float(e["value"]))
        report.append(dict(e, was=was, now=rec.get(e["field"])))
    # Rebuild from the end so earlier offsets stay valid.
    for rec in sorted({id(r): r for _, r in plan}.values(), key=lambda r: -r.start):
        text = text[:rec.start] + rec.text + text[rec.end:]
    Path(dst).write_text(text)
    return report
