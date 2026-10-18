# -------------------------------------------------------------
# file: gridpack/mcp/events.py
# -------------------------------------------------------------
# Dynamic events for a derived case.
# -------------------------------------------------------------

"""Faults and trips, written under ``Dynamic_simulation/Events``.

GridPACK applies an event only at a time step it lands on exactly, and
matches lines and generators by their RAW ids, so both are checked here.
"""

import xml.etree.ElementTree as ET

from ..rawedit import EditError, _id, _tokens, raw_blocks

__all__ = ["EVENT_OPS", "append_events", "plan_events", "raw_elements",
           "set_observations", "watch_all_generators"]

# op -> (required keys, optional keys)
EVENT_OPS = {
    "add_fault": (("from", "to", "start", "duration"), ("ckt",)),
    "add_bus_fault": (("bus", "start", "duration"), ("g_pu", "b_pu")),
    "trip_branch": (("from", "to", "time"), ("ckt", "in_service")),
    "trip_generator": (("bus", "time"), ("id", "in_service")),
}
_B_FAULT = 99999.0            # GridPACK's own default: a solid short


def raw_elements(path):
    """``(buses, branches, generators)`` of a RAW file: bus numbers,
    ``(from, to, ckt)`` including transformers, and ``(bus, id)``."""
    blocks = raw_blocks(path)
    buses = {int(_tokens(l)[0][0]) for _, l in blocks.get("bus", [])}
    gens = {(int(f[0]), _id(f[1])) for f in
            (_tokens(l)[0] for _, l in blocks.get("generator", []))}
    branches = set()
    for _, line in blocks.get("branch", []):
        f = _tokens(line)[0]
        branches.add((abs(int(f[0])), abs(int(f[1])), _id(f[2])))
    for _, line in blocks.get("transformer", []):
        # A record's first line: I, J, K, 'CKT', ...; winding lines differ.
        f = _tokens(line)[0]
        if len(f) > 3 and f[3].strip().startswith("'"):
            try:
                branches.add((abs(int(f[0])), abs(int(f[1])), _id(f[3])))
            except ValueError:
                pass
    return buses, branches, gens


def plan_events(edits, raw, time_step, end_time):
    """``[(edit, element)]`` for event edits, all checked against the RAW
    and the simulation's time step and length."""
    buses, branches, gens = raw_elements(raw)
    out = []
    for n, edit in enumerate(edits, 1):
        op = edit["op"]
        required, optional = EVENT_OPS[op]
        where = "edit %d (%s)" % (n, op)
        missing = [k for k in required if k not in edit]
        extra = [k for k in edit if k not in required + optional + ("op",)]
        if missing or extra:
            raise EditError("%s: %s" % (where, "; ".join(
                (["missing %s" % ", ".join(missing)] if missing else [])
                + (["unexpected %s" % ", ".join(extra)] if extra else []))))
        for k in ("bus", "from", "to"):
            if k in edit and (isinstance(edit[k], bool)
                              or not isinstance(edit[k], int)):
                raise EditError("%s: %s must be a bus number" % (where, k))
        for k in ("start", "duration", "time", "g_pu", "b_pu"):
            if k in edit and (isinstance(edit[k], bool)
                              or not isinstance(edit[k], (int, float))):
                raise EditError("%s: %s must be a number" % (where, k))
        if "in_service" in edit and not isinstance(edit["in_service"], bool):
            raise EditError("%s: in_service must be true or false" % where)
        if "bus" in edit and edit["bus"] not in buses:
            raise EditError("%s: no bus %d in the RAW" % (where, edit["bus"]))
        if "from" in edit:
            i, j, ckt = _branch(edit, branches, where)
        if op == "trip_generator":
            gid = _gen(edit, gens, where)
        times = {k: edit[k] for k in ("start", "time") if k in edit}
        if "duration" in edit:
            if edit["duration"] <= 0:
                raise EditError("%s: duration must be positive" % where)
            times["start + duration"] = edit["start"] + edit["duration"]
        for k, t in times.items():
            if t < 0:
                raise EditError("%s: %s must not be negative" % (where, k))
            if t > end_time + 1e-9:
                raise EditError("%s: %s = %g is past simulationTime %g"
                                % (where, k, t, end_time))
            if abs(t / time_step - round(t / time_step)) > 1e-6:
                raise EditError("%s: %s = %g is not a multiple of the time "
                                "step %g" % (where, k, t, time_step))
        status = "1" if edit.get("in_service", False) else "0"
        if op == "add_fault":
            el = _element("faultEvent", beginFault=edit["start"],
                          endFault=edit["start"] + edit["duration"],
                          faultBranch="%d %d" % (edit["from"], edit["to"]),
                          id=ckt)
        elif op == "add_bus_fault":
            el = _element("BusFault", begin=edit["start"],
                          end=edit["start"] + edit["duration"],
                          bus=edit["bus"], yfault="%g %g" % (
                              edit.get("g_pu", 0.0), edit.get("b_pu", _B_FAULT)))
        elif op == "trip_branch":
            el = _element("LineStatus", time=edit["time"],
                          line="%d %d" % (i, j), id=ckt, status=status)
        else:
            el = _element("GenStatus", time=edit["time"], bus=edit["bus"],
                          id=gid, status=status)
        out.append((edit, el))
    return out


def append_events(xml, elements, replace=False):
    """Append ``elements`` under Dynamic_simulation/Events of ``xml``;
    with ``replace`` the existing events go first."""
    tree = ET.parse(xml)
    ds = tree.getroot().find("Dynamic_simulation")
    if ds is None:
        raise EditError("%s has no Dynamic_simulation block" % xml)
    events = ds.find("Events")
    if events is None:
        events = ET.SubElement(ds, "Events")
    if replace:
        for old in list(events):
            events.remove(old)
    events.extend(elements)
    ET.indent(tree.getroot(), "  ")
    tree.write(xml, encoding="utf-8", xml_declaration=True)


def watch_all_generators(xml, raw):
    """Point Dynamic_simulation/generatorWatch of ``xml`` at every
    generator in ``raw``; return how many."""
    _, _, gens = raw_elements(raw)
    tree = ET.parse(xml)
    ds = tree.getroot().find("Dynamic_simulation")
    if ds is None:
        raise EditError("%s has no Dynamic_simulation block" % xml)
    watch = ds.find("generatorWatch")
    if watch is None:
        watch = ET.SubElement(ds, "generatorWatch")
    for old in list(watch):
        watch.remove(old)
    for bus, gid in sorted(gens):
        g = ET.SubElement(watch, "generator")
        ET.SubElement(g, "busID").text = str(bus)
        ET.SubElement(g, "generatorID").text = gid
    if ds.find("generatorWatchFileName") is None:
        ET.SubElement(ds, "generatorWatchFileName").text = "gen_watch.csv"
    ET.indent(tree.getroot(), "  ")
    tree.write(xml, encoding="utf-8", xml_declaration=True)
    return len(gens)


def set_observations(xml, raw, generators=None, buses=None,
                     bus_frequency=False):
    """Replace Dynamic_simulation/observations of ``xml`` (what HADREC
    records): ``generators`` and ``buses`` are "all" or bus-number lists;
    ``bus_frequency`` adds a frequency observation per listed bus.
    Returns the counts by type."""
    all_buses, _, gens = raw_elements(raw)
    tree = ET.parse(xml)
    ds = tree.getroot().find("Dynamic_simulation")
    if ds is None:
        raise EditError("%s has no Dynamic_simulation block" % xml)
    obs = ds.find("observations")
    if obs is None:
        obs = ET.SubElement(ds, "observations")
    for old in list(obs):
        obs.remove(old)
    counts = {}

    def add(kind, bus, gid=None):
        o = ET.SubElement(obs, "observation")
        ET.SubElement(o, "type").text = kind
        ET.SubElement(o, "busID").text = str(bus)
        if gid is not None:
            ET.SubElement(o, "generatorID").text = gid
        counts[kind] = counts.get(kind, 0) + 1

    if generators == "all":
        for bus, gid in sorted(gens):
            add("generator", bus, gid)
    elif generators:
        for bus in generators:
            mine = sorted(g for b, g in gens if b == bus)
            if not mine:
                raise EditError("no generator at bus %d" % bus)
            for gid in mine:
                add("generator", bus, gid)
    chosen = sorted(all_buses) if buses == "all" else list(buses or [])
    for bus in chosen:
        if bus not in all_buses:
            raise EditError("no bus %d in the RAW" % bus)
        add("bus", bus)
    if bus_frequency:
        for bus in chosen:
            add("busfrequency", bus)
    ET.indent(tree.getroot(), "  ")
    tree.write(xml, encoding="utf-8", xml_declaration=True)
    return counts


def _branch(edit, branches, where):
    """``(from, to, ckt)`` as the RAW writes the branch."""
    i, j = edit["from"], edit["to"]
    found = {(a, b, c) for a, b, c in branches if {a, b} == {i, j}}
    ckts = sorted(c for _, _, c in found)
    if not found:
        raise EditError("%s: no branch %d-%d in the RAW" % (where, i, j))
    if "ckt" in edit:
        ckt = _id(str(edit["ckt"]))
        if ckt not in ckts:
            raise EditError("%s: branch %d-%d has circuit %s, not %s" % (
                where, i, j, ", ".join(ckts), ckt))
    elif len(found) > 1:
        raise EditError("%s: branch %d-%d has circuits %s; give ckt" % (
            where, i, j, ", ".join(ckts)))
    else:
        ckt = ckts[0]
    return next(b for b in found if b[2] == ckt)


def _gen(edit, gens, where):
    found = {g for b, g in gens if b == edit["bus"]}
    if not found:
        raise EditError("%s: no generator at bus %d" % (where, edit["bus"]))
    if "id" in edit:
        gid = _id(str(edit["id"]))
        if gid not in found:
            raise EditError("%s: bus %d has generator id %s, not %s" % (
                where, edit["bus"], ", ".join(sorted(found)), gid))
        return gid
    if len(found) > 1:
        raise EditError("%s: bus %d has generators %s; give id" % (
            where, edit["bus"], ", ".join(sorted(found))))
    return found.pop()


def _element(tag, **children):
    el = ET.Element(tag)
    for name, value in children.items():
        ET.SubElement(el, name).text = "%g" % value \
            if isinstance(value, float) else str(value)
    return el
