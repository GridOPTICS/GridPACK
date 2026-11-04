# -------------------------------------------------------------
# file: gridpack/rawedit.py
# -------------------------------------------------------------
# Typed edits to the records of a RAW file.
# -------------------------------------------------------------

"""Change load, generator, shunt and branch records in a RAW file.

Only the named fields of the named records change; every other byte of
the file is kept.  Supported for RAW versions 23, 33, 34, 35 and 36.
"""

import json
import re
from pathlib import Path

__all__ = ["OPS", "edit_raw", "raw_blocks", "raw_generators", "raw_version"]

_END = re.compile(r"\s*(0|Q)(\s|/|,|$)")
_BEGIN = re.compile(r"BEGIN (.+?) DATA", re.I)
_LABELS = {"BUS": "bus", "LOAD": "load", "FIXED SHUNT": "fixed_shunt",
           "GENERATOR": "generator", "BRANCH": "branch", "LINE": "branch",
           "TRANSFORMER": "transformer"}
_ORDER = {23: ("bus", "generator", "branch", "transformer_adjustment"),
          33: ("bus", "load", "fixed_shunt", "generator", "branch",
               "transformer")}

# Field names per record kind; the entry for the highest version <= the
# file's applies.  v23 keeps loads and shunts in the bus record.
_FIELDS = {
    "bus": {23: "I IDE PL QL GL BL IA VM VA NAME BASKV ZONE",
            33: "I NAME BASKV IDE AREA ZONE OWNER VM VA NVHI NVLO EVHI EVLO"},
    "load": {33: "I ID STAT AREA ZONE PL QL IP IQ YP YQ OWNER SCALE INTRPT "
                 "DGENP DGENQ DGENF"},
    "fixed_shunt": {33: "I ID STAT GL BL"},
    "generator": {23: "I ID PG QG QT QB VS IREG MBASE ZR ZX RT XT GTAP STAT "
                      "RMPCT PT PB",
                  35: "I ID PG QG QT QB VS IREG NREG MBASE ZR ZX RT XT GTAP "
                      "STAT RMPCT PT PB"},
    "branch": {23: "I J CKT R X B RATEA RATEB RATEC RATIO ANGLE GI BI GJ BJ "
                   "STAT",
               33: "I J CKT R X B RATEA RATEB RATEC GI BI GJ BJ STAT MET LEN",
               34: "I J CKT R X B NAME RATE1 RATE2 RATE3 RATE4 RATE5 RATE6 "
                   "RATE7 RATE8 RATE9 RATE10 RATE11 RATE12 GI BI GJ BJ STAT "
                   "MET LEN"},
}
_INT_FIELDS = {"I", "J", "IDE", "IA", "ZONE", "AREA", "OWNER", "STAT",
               "IREG", "NREG", "MET", "INTRPT", "SCALE"}

# op -> (record kind, required keys, optional keys)
OPS = {
    "set_load_active_power": ("load", ("bus", "mw"), ("id",)),
    "set_load_reactive_power": ("load", ("bus", "mvar"), ("id",)),
    "scale_load": ("load", ("bus", "factor"), ("id",)),
    "set_load_in_service": ("load", ("bus", "in_service"), ("id",)),
    "set_generator_active_power": ("generator", ("bus", "mw"), ("id",)),
    "set_generator_voltage_setpoint": ("generator", ("bus", "vm_pu"), ("id",)),
    "set_generator_in_service": ("generator", ("bus", "in_service"), ("id",)),
    "set_shunt": ("fixed_shunt", ("bus",), ("id", "mvar", "mw")),
    "set_shunt_in_service": ("fixed_shunt", ("bus", "in_service"), ("id",)),
    "set_branch_in_service": ("branch", ("from", "to", "ckt", "in_service"),
                              ()),
    "scale_all_loads": ("load", ("factor",), ()),
    "scale_all_generation": ("generator", ("factor",), ()),
    "set_param": (None, ("kind", "field", "value"),
                  ("bus", "id", "from", "to", "ckt")),
}
_KINDS = ("bus", "load", "generator", "fixed_shunt", "branch")


class EditError(ValueError):
    """An edit that cannot be applied; nothing has been written."""


# ---- reading -----------------------------------------------------------

def _lines(path):
    # newline="" keeps the file's own line endings for the rewrite.
    with open(path, newline="", errors="replace") as fh:
        return fh.read().splitlines(keepends=True)


def raw_version(path) -> int:
    """The RAW revision: the header's third field, else 23."""
    for line in _lines(path):
        if not line.lstrip().startswith("@!"):
            fields = line.split("/")[0].split(",")
            try:
                return int(float(fields[2])) if len(fields) >= 3 else 23
            except ValueError:
                return 23
    return 23


def raw_blocks(path):
    """``{kind: [(line number, line)]}`` for the data blocks of a RAW file.

    Blocks are named by the ``BEGIN ... DATA`` labels when the file has
    them, else by position for the version.
    """
    lines = _lines(path)
    version = raw_version(path)
    order = _ORDER[23 if version < 33 else 33]
    data = [(i, l) for i, l in enumerate(lines) if not l.lstrip().startswith("@!")]
    blocks, cur, before, n = {}, [], None, 0
    for i, line in data[3:]:
        if not line.strip():
            continue
        if _END.match(line):
            m = _BEGIN.search(line)
            after = _LABELS.get(m.group(1).upper()) if m else None
            if before:
                name = before
            elif after == "bus" and n == 0:
                name = "system"          # v34+ system-wide block first
            else:
                name = order[n] if n < len(order) else "block%d" % n
            blocks.setdefault(name, []).extend(cur)
            cur, before, n = [], after, n + 1
        else:
            cur.append((i, line))
    return blocks


def raw_generators(path):
    """``[(bus, id)]`` for every generator record in a RAW file."""
    out = []
    for _, line in raw_blocks(path).get("generator", []):
        tok, _, _ = _tokens(line)
        try:
            out.append((int(tok[0]), _id(tok[1])))
        except (IndexError, ValueError):
            raise ValueError("unreadable generator record in %s: %r"
                             % (Path(path).name, line.strip()))
    return out


def _tokens(line):
    """``(fields, delimiter, tail)``: the record's fields outside quotes,
    the delimiter used, and the comment plus line ending."""
    body = line.rstrip("\r\n")
    end = line[len(body):]
    out, cur, quoted = [], [], False
    for i, ch in enumerate(body):
        if ch == "'":
            quoted = not quoted
        elif ch == "/" and not quoted:
            out.append("".join(cur))
            return _resplit(out, body[i:] + end)
        elif ch == "," and not quoted:
            out.append("".join(cur))
            cur = []
            continue
        cur.append(ch)
    out.append("".join(cur))
    return _resplit(out, end)


def _resplit(fields, tail):
    if len(fields) == 1 and len(fields[0].split()) > 1:
        return fields[0].split(), " ", tail      # space-delimited v23
    return fields, ",", tail


def _id(token):
    return token.strip().strip("'\"").strip() or "1"


def _fields(kind, version):
    table = _FIELDS[kind]
    v = max(k for k in table if k <= version)
    return table[v].split()


# ---- editing -----------------------------------------------------------

def _fmt(value, old):
    """``value`` written like ``old``: its decimals (more if needed to
    hold the value) and its width, so columns stay aligned."""
    if isinstance(value, bool):
        value = int(value)
    if isinstance(value, int):
        text = "%d" % value
    elif value != 0 and abs(value) < 1e-5:
        text = "%.6e" % value
    else:
        needed = len(("%.6f" % value).rstrip("0").split(".")[1])
        had = len(old.strip().split(".")[1]) if "." in old and \
            "E" not in old.upper() else 0
        text = "%.*f" % (max(needed, had, 1), value)
    return text.rjust(len(old)) if len(text) < len(old) else text


def _number(token):
    text = token.strip()
    try:
        return int(text)
    except ValueError:
        return float(text)


class _Record:
    """One data line, editable by field name."""

    def __init__(self, lineno, line, names):
        self.lineno = lineno
        self.fields, self.delim, self.tail = _tokens(line)
        self.names = names

    def index(self, name):
        if name not in self.names:
            raise EditError("field %s is not in a %s record" % (name, "/".join(
                self.names[:3])))
        i = self.names.index(name)
        if i >= len(self.fields):
            raise EditError("record has no %s field: %r" % (name, self.text()))
        return i

    def get(self, name):
        return _number(self.fields[self.index(name)])

    def set(self, name, value):
        i = self.index(name)
        if name in _INT_FIELDS:
            value = int(value)
        self.fields[i] = _fmt(value, self.fields[i])

    def text(self):
        return self.delim.join(self.fields) + self.tail


class _Raw:
    """A RAW file's records, indexed for lookup."""

    def __init__(self, path):
        self.path = Path(path)
        self.lines = _lines(path)
        self.version = raw_version(path)
        self.blocks = raw_blocks(path)
        self.records = {}
        for kind in _KINDS:
            if kind in ("load", "fixed_shunt") and self.version < 33:
                continue
            names = _fields(kind, self.version)
            for lineno, line in self.blocks.get(kind, []):
                rec = _Record(lineno, line, names)
                self.records[(kind, self._key(kind, rec))] = rec

    @staticmethod
    def _key(kind, rec):
        f = rec.fields
        if kind == "bus":
            return int(f[0])
        if kind == "branch":
            return (abs(int(f[0])), abs(int(f[1])), _id(f[2]))
        return (int(f[0]), _id(f[1]))

    def find(self, kind, key, what):
        if kind == "branch":
            i, j, ckt = key
            rec = self.records.get((kind, (i, j, ckt))) or \
                self.records.get((kind, (j, i, ckt)))
        else:
            rec = self.records.get((kind, key))
        if rec is None:
            hint = ("; transformers are in their own block and cannot be "
                    "switched yet" if kind == "branch" else "")
            raise EditError("no %s %s in %s%s" % (kind, what, self.path.name,
                                                 hint))
        return rec

    def write(self, dst):
        for rec in self.records.values():
            self.lines[rec.lineno] = rec.text()
        with open(dst, "w", newline="") as fh:
            fh.writelines(self.lines)


def _check(edit, n):
    if not isinstance(edit, dict) or "op" not in edit:
        raise EditError("edit %d must be an object with an \"op\"" % n)
    op = edit["op"]
    if op not in OPS:
        raise EditError("edit %d: unknown op %r; one of %s"
                        % (n, op, ", ".join(OPS)))
    kind, required, optional = OPS[op]
    missing = [k for k in required if k not in edit]
    extra = [k for k in edit if k not in required + optional + ("op",)]
    if missing or extra:
        raise EditError("edit %d (%s): %s" % (n, op, "; ".join(
            (["missing %s" % ", ".join(missing)] if missing else [])
            + (["unexpected %s" % ", ".join(extra)] if extra else []))))
    for k in ("bus", "from", "to"):
        if k in edit and (isinstance(edit[k], bool)
                          or not isinstance(edit[k], int)):
            raise EditError("edit %d (%s): %s must be a bus number" % (n, op, k))
    for k in ("mw", "mvar", "factor", "vm_pu"):
        if k in edit and (isinstance(edit[k], bool)
                          or not isinstance(edit[k], (int, float))):
            raise EditError("edit %d (%s): %s must be a number" % (n, op, k))
    if "in_service" in edit and not isinstance(edit["in_service"], bool):
        raise EditError("edit %d (%s): in_service must be true or false"
                        % (n, op))
    if op == "set_param" and edit["kind"] not in _KINDS:
        raise EditError("edit %d: kind must be one of %s" % (n, ", ".join(_KINDS)))
    if op == "set_param" and (isinstance(edit["value"], bool) or
                              not isinstance(edit["value"], (int, float))):
        raise EditError("edit %d (set_param): value must be a number" % n)


def _plan(raw, edit, n):
    """``[(record, fields, fn)]`` for one validated edit: ``fn(record)``
    gives the new field values when the edit is applied, so earlier edits
    to the same record are seen."""
    op = edit["op"]
    kind = OPS[op][0] or edit["kind"]
    v23 = raw.version < 33
    if op in ("scale_all_loads", "scale_all_generation"):
        f = edit["factor"]
        if op == "scale_all_loads":
            recs = [r for (k, _), r in raw.records.items()
                    if k == ("bus" if v23 else "load")]
            return [(r, ["PL", "QL"], lambda r, f=f: {
                "PL": r.get("PL") * f, "QL": r.get("QL") * f}) for r in recs]
        recs = [r for (k, _), r in raw.records.items() if k == "generator"]
        return [(r, ["PG"], lambda r, f=f: {"PG": r.get("PG") * f})
                for r in recs]
    rec, new = _plan_one(raw, edit, n, op, kind, v23)
    if callable(new):
        return [(rec, new.names, new)]
    return [(rec, list(new), lambda r, v=new: v)]


def _plan_one(raw, edit, n, op, kind, v23):
    ident = _id(str(edit.get("id", "1")))
    if kind == "branch":
        key = (edit["from"], edit["to"], _id(str(edit["ckt"])))
        what = "%d-%d circuit %s" % key
    elif kind == "bus":
        key, what = edit["bus"], "%d" % edit["bus"]
    else:
        key, what = (edit["bus"], ident), "%d id %s" % (edit["bus"], ident)
    if v23 and kind in ("load", "fixed_shunt"):
        # One load and one shunt per bus, inside the bus record.
        if ident != "1":
            raise EditError("edit %d: a v23 RAW has one %s per bus; id must "
                            "be 1" % (n, kind))
        if op.endswith("_in_service"):
            raise EditError("edit %d: a v23 RAW has no %s status; set its "
                            "power to 0 instead" % (n, kind))
        rec = raw.find("bus", edit["bus"], "%d" % edit["bus"])
    else:
        rec = raw.find(kind, key, what)
    if op == "set_load_active_power":
        return rec, {"PL": edit["mw"]}
    if op == "set_load_reactive_power":
        return rec, {"QL": edit["mvar"]}
    if op == "scale_load":
        f = edit["factor"]
        fn = lambda r, f=f: {"PL": r.get("PL") * f, "QL": r.get("QL") * f}
        fn.names = ["PL", "QL"]
        return rec, fn
    if op == "set_generator_active_power":
        return rec, {"PG": edit["mw"]}
    if op == "set_generator_voltage_setpoint":
        return rec, {"VS": edit["vm_pu"]}
    if op == "set_shunt":
        new = {}
        if "mvar" in edit:
            new["BL"] = edit["mvar"]
        if "mw" in edit:
            new["GL"] = edit["mw"]
        if not new:
            raise EditError("edit %d (set_shunt): give mvar and/or mw" % n)
        return rec, new
    if op.endswith("_in_service"):
        return rec, {"STAT": 1 if edit["in_service"] else 0}
    return rec, {edit["field"]: edit["value"]}           # set_param


def edit_raw(src, dst, edits):
    """Apply ``edits`` (a list of op dicts, see :data:`OPS`) to the RAW
    file ``src`` and write ``dst``; return what each edit changed.

    Every edit is checked before anything is written.
    """
    if not isinstance(edits, list) or not edits:
        raise EditError("edits must be a non-empty list")
    raw = _Raw(src)
    plan = []
    for n, edit in enumerate(edits, 1):
        _check(edit, n)
        targets = _plan(raw, edit, n)
        if not targets:
            raise EditError("edit %d (%s): no records to change" % (n, edit["op"]))
        for rec, names, _ in targets:
            for name in names:
                rec.index(name)               # field exists in this version
        plan.append((edit, targets))
    report = []
    for edit, targets in plan:
        was = _values(targets)
        for rec, _, fn in targets:
            for name, value in fn(rec).items():
                rec.set(name, value)
        entry = dict(edit, was=was, now=_values(targets))
        if len(targets) > 1:
            entry["records"] = len(targets)
        report.append(entry)
    raw.write(dst)
    return report


def _values(targets):
    """Field values of one record, or their totals over many."""
    if len(targets) == 1:
        rec, names, _ = targets[0]
        return {n: rec.get(n) for n in names}
    return {n: round(sum(rec.get(n) for rec, _, _ in targets), 6)
            for n in targets[0][1]}
