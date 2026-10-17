# -------------------------------------------------------------
# file: gridpack/_observations.py
# -------------------------------------------------------------
# Rebuild the global observation set from per-rank slices.
#
# DSFullApp::setObservations sums its per-observation "active" flags with
# a count of p_comm.size() instead of the observation count, so past the
# first nranks entries each rank lists, and returns values for, only the
# elements it owns.  This merges the slices back into XML order.
# -------------------------------------------------------------

from __future__ import annotations

import xml.etree.ElementTree as ET
from typing import Dict, List, Optional, Sequence, Tuple

# Field order of getObservationLists(_withBusFreq).
_LISTS = ("gen_buses", "gen_ids", "load_buses", "load_ids", "bus_ids",
          "bus_freq_ids")


def _norm(tag) -> str:
    return str(tag).replace("'", "").replace('"', "").strip()


def _xml_order(input_file: str) -> Optional[Dict[str, list]]:
    """Observation keys per kind in document order, as setObservations reads."""
    try:
        with open(input_file, "rb") as fh:
            root = ET.fromstring(fh.read().lstrip())
    except (OSError, ET.ParseError):
        return None
    out = {"gen": [], "bus": [], "load": [], "freq": []}
    for obs in root.iter("observation"):
        kind = (obs.findtext("type") or "").strip().lower()
        bus = (obs.findtext("busID") or "").strip()
        if not bus:
            continue
        try:
            bus = int(bus)
        except ValueError:
            continue
        if kind == "generator" and obs.findtext("generatorID") is not None:
            out["gen"].append((bus, _norm(obs.findtext("generatorID"))))
        elif kind == "bus":
            out["bus"].append(bus)
        elif kind == "busfrequency":
            out["freq"].append(bus)
        elif kind == "load" and obs.findtext("loadID") is not None:
            out["load"].append((bus, _norm(obs.findtext("loadID"))))
    return out


class ObservationMerger:
    """Per-rank observation lists and values -> the global set.

    Collective: construct and :meth:`merge` on every rank.  A no-op at one
    rank or without mpi4py.
    """

    def __init__(self, mpi_comm, input_file: str, lists: Sequence) -> None:
        lists = list(lists) + [[]] * (len(_LISTS) - len(lists))
        local = {k: list(v) for k, v in zip(_LISTS, lists)}
        self.active = mpi_comm is not None and mpi_comm.Get_size() > 1
        self._comm = mpi_comm
        if not self.active:
            self.lists = local
            return

        per_rank = mpi_comm.allgather(local)
        keys = [self._keys(r) for r in per_rank]
        order = _xml_order(input_file)

        # Global key order per kind: XML order, kept only where some rank
        # found the element; first-seen order if the XML cannot be read.
        glob: Dict[str, list] = {}
        for kind in ("gen", "bus", "load", "freq"):
            seen = {k for r in keys for k in r[kind]}
            if order is not None:
                g = [k for k in order[kind] if k in seen]
                # Anything the XML pass missed still gets a column.
                g += [k for r in keys for k in r[kind] if k not in g]
            else:
                g = []
                for r in keys:
                    g += [k for k in r[kind] if k not in g]
            glob[kind] = list(dict.fromkeys(g))

        # The ids as the C++ spells them, for the merged lists.
        spelled: Dict[Tuple, Tuple] = {}
        for r in per_rank:
            for b, i in zip(r["gen_buses"], r["gen_ids"]):
                spelled.setdefault(("g", b, _norm(i)), (b, i))
            for b, i in zip(r["load_buses"], r["load_ids"]):
                spelled.setdefault(("l", b, _norm(i)), (b, i))
        self.lists = {
            "gen_buses": [spelled[("g",) + k][0] for k in glob["gen"]],
            "gen_ids": [spelled[("g",) + k][1] for k in glob["gen"]],
            "load_buses": [spelled[("l",) + k][0] for k in glob["load"]],
            "load_ids": [spelled[("l",) + k][1] for k in glob["load"]],
            "bus_ids": list(glob["bus"]),
            "bus_freq_ids": list(glob["freq"]),
        }

        self._n = {k: len(v) for k, v in glob.items()}
        # Per rank: (global offset, local offset) pairs over the flat vector.
        self._maps: List[List[Tuple[int, int]]] = []
        for r in keys:
            pos = {kind: {k: j for j, k in enumerate(glob[kind])}
                   for kind in glob}
            n_loc = {kind: len(r[kind]) for kind in r}
            pairs = []
            g_off = l_off = 0
            # rSpd, rAng, P, Q per generator; vMag, vAng per bus; then
            # fOnline per load and frequency per bus-frequency entry.
            for kind, reps in (("gen", 4), ("bus", 2), ("load", 1),
                               ("freq", 1)):
                for _ in range(reps):
                    for j, k in enumerate(r[kind]):
                        pairs.append((g_off + pos[kind][k], l_off + j))
                    g_off += self._n[kind]
                    l_off += n_loc[kind]
            self._maps.append(pairs)
        self._width = g_off

    @staticmethod
    def _keys(r: Dict[str, list]) -> Dict[str, list]:
        return {
            "gen": [(b, _norm(i)) for b, i in zip(r["gen_buses"], r["gen_ids"])],
            "bus": list(r["bus_ids"]),
            "load": [(b, _norm(i)) for b, i in zip(r["load_buses"], r["load_ids"])],
            "freq": list(r["bus_freq_ids"]),
        }

    def as_tuple(self, with_bus_freq: bool) -> tuple:
        names = _LISTS if with_bus_freq else _LISTS[:-1]
        return tuple(self.lists[k] for k in names)

    def merge(self, obs: Sequence[float]) -> list:
        """This step's global observation vector.  Collective."""
        if not self.active:
            return obs
        chunks = self._comm.allgather([float(v) for v in obs])
        out = [float("nan")] * self._width
        for pairs, chunk in zip(self._maps, chunks):
            for g, l in pairs:
                if l < len(chunk):
                    out[g] = chunk[l]
        return out
