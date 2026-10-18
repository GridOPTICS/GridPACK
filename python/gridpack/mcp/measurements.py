# -------------------------------------------------------------
# file: gridpack/mcp/measurements.py
# -------------------------------------------------------------
# Synthetic state-estimation measurements from a power flow solution.
# -------------------------------------------------------------

"""Meters, noise and planted bad data for a state estimation case.

Values are what GridPACK's estimator reads: voltages in p.u., angles in
degrees, powers in p.u. on the case base.  The noise standard deviation
goes into each meter's ``Deviation`` so the estimator weights it right.
"""

import json
import math
import random
import xml.etree.ElementTree as ET

__all__ = ["DEFAULT_STD", "DEFAULT_TYPES", "TYPES", "accuracy", "dump",
           "estimate_lines", "raw_sbase", "se_xml", "synthesize",
           "write_measurements"]

TYPES = ("VM", "VA", "PI", "QI", "PIJ", "QIJ", "PJI", "QJI")
DEFAULT_TYPES = ("VM", "PI", "QI", "PIJ", "QIJ")
# Typical meter accuracies: p.u. except VA in degrees.
DEFAULT_STD = {"VM": 0.004, "VA": 0.1, "PI": 0.01, "QI": 0.01,
               "PIJ": 0.008, "QIJ": 0.008, "PJI": 0.008, "QJI": 0.008}
_BUS = {"VM": "voltage", "VA": "angle", "PI": "pInjection", "QI": "qInjection"}
_BRANCH = {"PIJ": "pFrom", "QIJ": "qFrom", "PJI": "pTo", "QJI": "qTo"}
_PU = {"PI", "QI", "PIJ", "QIJ", "PJI", "QJI"}

_SE_XML = """<?xml version="1.0" encoding="utf-8"?>
<Configuration>
  <State_estimation>
    <%(tag)s>%(raw)s</%(tag)s>
    <measurementList>%(meas)s</measurementList>
    <tolerance>1.0e-3</tolerance>
    <maxIteration>20</maxIteration>
    <badDataThreshold>3.0</badDataThreshold>
    <useVoltageConstraints>false</useVoltageConstraints>
    <enforceVoltageLimits>false</enforceVoltageLimits>
    <minVoltage>0.9</minVoltage>
    <maxVoltage>1.1</maxVoltage>
    <voltageConstraintDeviation>0.0025</voltageConstraintDeviation>
    <LinearSolver>
      <PETScOptions>
        -ksp_type richardson
        -pc_type lu
        -pc_factor_mat_solver_type superlu_dist
        -ksp_max_it 1
      </PETScOptions>
    </LinearSolver>
  </State_estimation>
</Configuration>
"""


def raw_sbase(path):
    """The system base MVA from the RAW header (second field), else 100."""
    with open(path, errors="replace") as fh:
        for line in fh:
            if line.lstrip().startswith("@!"):
                continue
            body = line.split("/")[0]
            fields = body.split(",") if "," in body else body.split()
            try:
                return float(fields[1])
            except (IndexError, ValueError):
                return 100.0
    return 100.0


def se_xml(raw, version, meas):
    """A State_estimation input around ``raw`` and the meter file."""
    tag = "networkConfiguration" if version == 23 \
        else "networkConfiguration_v%d" % version
    return _SE_XML % {"tag": tag, "raw": raw, "meas": meas}


def _key(m):
    if "bus" in m:
        return (m["type"], m["bus"])
    return (m["type"], m["from"], m["to"], m["ckt"])


def _truth(pf, sbase, types):
    out = []
    for t in types:
        if t in _BUS:
            for b in pf["buses"]:
                v = b[_BUS[t]]
                out.append({"type": t, "bus": b["busId"],
                            "truth": v / sbase if t in _PU else v})
        else:
            for br in pf["branches"]:
                out.append({"type": t, "from": br["fromBus"], "to": br["toBus"],
                            "ckt": br["circuitId"].strip(),
                            "truth": br[_BRANCH[t]] / sbase})
    return out


def synthesize(pf, sbase, types=None, std=None, noise_scale=1.0, coverage=1.0,
               bad_data=(), random_bad=0, bad_sigma=20.0, seed=0):
    """The meter list for a measurement file, plus a summary.

    Each meter has its truth, the value written (truth + noise [+ planted
    error]) and its deviation (the noise std).  ``coverage`` keeps that
    fraction of meters at random; planted meters are always kept.  With
    ``noise_scale`` 0 the values are exact and the deviations nominal.
    """
    types = list(types or DEFAULT_TYPES)
    unknown = [t for t in types if t not in TYPES]
    if unknown:
        raise ValueError("unknown measurement type %s; types are %s"
                         % (", ".join(unknown), ", ".join(TYPES)))
    std = dict(DEFAULT_STD, **(std or {}))
    if not 0 < coverage <= 1:
        raise ValueError("coverage must be in (0, 1]")
    rng = random.Random(seed)
    meters = _truth(pf, sbase, types)
    by_key = {_key(m): m for m in meters}
    planted = []
    for n, bad in enumerate(bad_data or (), 1):
        key = _bad_key(bad, n)
        if key not in by_key:
            raise ValueError("bad_data %d: no %s meter %s in the synthesized "
                             "set (type not chosen, or no such element)"
                             % (n, key[0], key[1:]))
        planted.append((by_key[key], float(bad["error"])))
    if coverage < 1:
        keep = set(id(m) for m, _ in planted)
        pool = [m for m in meters if id(m) not in keep]
        n_keep = max(0, round(coverage * len(meters)) - len(keep))
        keep |= set(id(m) for m in rng.sample(pool, min(n_keep, len(pool))))
        meters = [m for m in meters if id(m) in keep]
    if random_bad:
        if random_bad > len(meters):
            raise ValueError("random_bad %d exceeds the %d meters"
                             % (random_bad, len(meters)))
        taken = set(id(m) for m, _ in planted)
        pool = [m for m in meters if id(m) not in taken]
        for m in rng.sample(pool, min(random_bad, len(pool))):
            planted.append((m, rng.choice((-1, 1)) * bad_sigma * std[m["type"]]))
    errors = {id(m): e for m, e in planted}
    for m in meters:
        s = std[m["type"]] * noise_scale
        # Deviation is also the estimator's weight: exact values keep the
        # nominal one, or every rounding error reads as bad data.
        m["deviation"] = s if s > 0 else std[m["type"]]
        m["value"] = m["truth"] + (rng.gauss(0.0, s) if s > 0 else 0.0) \
            + errors.get(id(m), 0.0)
        if id(m) in errors:
            m["planted_error"] = errors[id(m)]
    counts = {}
    for m in meters:
        counts[m["type"]] = counts.get(m["type"], 0) + 1
    return meters, {"meters": len(meters), "by_type": counts,
                    "planted": [dict({k: v for k, v in m.items()
                                      if k not in ("value", "deviation", "truth")})
                                for m, _ in planted],
                    "seed": seed, "noise_scale": noise_scale,
                    "coverage": coverage, "std": {t: std[t] for t in types}}


def _bad_key(bad, n):
    if not isinstance(bad, dict) or "type" not in bad or "error" not in bad:
        raise ValueError("bad_data %d must have type and error" % n)
    t = bad["type"]
    if t in _BUS:
        if "bus" not in bad:
            raise ValueError("bad_data %d (%s) needs bus" % (n, t))
        return (t, bad["bus"])
    if t in _BRANCH:
        for k in ("from", "to", "ckt"):
            if k not in bad:
                raise ValueError("bad_data %d (%s) needs from, to and ckt" % (n, t))
        return (t, bad["from"], bad["to"], str(bad["ckt"]).strip())
    raise ValueError("bad_data %d: unknown type %r" % (n, t))


def write_measurements(path, meters):
    """Write the ``<Measurements>`` file GridPACK's estimator reads."""
    root = ET.Element("Measurements")
    for m in meters:
        el = ET.SubElement(root, "Measurement")
        ET.SubElement(el, "Type").text = m["type"]
        if "bus" in m:
            ET.SubElement(el, "Bus").text = str(m["bus"])
        else:
            ET.SubElement(el, "FromBus").text = str(m["from"])
            ET.SubElement(el, "ToBus").text = str(m["to"])
            ET.SubElement(el, "CKT").text = m["ckt"]
        ET.SubElement(el, "Value").text = "%.6f" % m["value"]
        ET.SubElement(el, "Deviation").text = "%.6g" % m["deviation"]
    ET.indent(root, "  ")
    ET.ElementTree(root).write(path, encoding="utf-8", xml_declaration=True)


def accuracy(se, truth_buses, planted, flagged):
    """Estimate against the power flow truth, and planted against flagged
    bad data."""
    truth = {b["busId"]: b for b in truth_buses}
    est = [b for b in se.get("buses", []) if b["busId"] in truth]
    out = {"buses": len(est)}
    for key, unit in (("voltage", "pu"), ("angle", "deg")):
        errs = [(b[key] - truth[b["busId"]][key], b["busId"]) for b in est]
        if errs:
            worst = max(errs, key=lambda x: abs(x[0]))
            out["%s_%s" % (key, unit)] = {
                "rms_error": math.sqrt(sum(e * e for e, _ in errs) / len(errs)),
                "max_error": worst[0], "at_bus": worst[1]}
    want = {_key(m) for m in planted}
    got = {_key(_norm(f)) for f in flagged}
    out["bad_data"] = {"planted": len(want), "flagged": len(got),
                       "caught": sorted(" ".join(str(x) for x in k) for k in want & got),
                       "missed": sorted(" ".join(str(x) for x in k) for k in want - got),
                       "false_alarms": sorted(" ".join(str(x) for x in k) for k in got - want)}
    return out


def _norm(f):
    if "busId" in f:
        return {"type": f["type"], "bus": f["busId"]}
    return {"type": f["type"], "from": f["fromBus"], "to": f["toBus"],
            "ckt": str(f.get("ckt", f.get("circuitId", ""))).strip()}


def estimate_lines(acc):
    """Report lines for :func:`accuracy`."""
    lines = []
    if "voltage_pu" in acc:
        v, a = acc["voltage_pu"], acc["angle_deg"]
        lines.append("Against the power flow truth: RMS voltage error %.5f p.u. "
                     "(largest %+.5f at bus %d), RMS angle error %.4f deg "
                     "(largest %+.4f at bus %d)."
                     % (v["rms_error"], v["max_error"], v["at_bus"],
                        a["rms_error"], a["max_error"], a["at_bus"]))
    bd = acc["bad_data"]
    if bd["planted"] or bd["flagged"]:
        lines.append("Bad data: %d planted, %d flagged; caught %s; missed %s; "
                     "false alarms %s." % (
                         bd["planted"], bd["flagged"],
                         ", ".join(bd["caught"]) or "none",
                         ", ".join(bd["missed"]) or "none",
                         ", ".join(bd["false_alarms"]) or "none"))
    return lines


def dump(path, pf, summary):
    """Save the truth and the synthesis settings beside the case."""
    path.write_text(json.dumps({"truth_buses": pf["buses"], **summary}, indent=1))
