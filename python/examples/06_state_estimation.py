#!/usr/bin/env python
# -------------------------------------------------------------
# file: python/examples/06_state_estimation.py
# -------------------------------------------------------------

"""Estimate the IEEE 14 state from noisy meters, then corrupt one.

    python 06_state_estimation.py
    mpiexec -np 2 python 06_state_estimation.py

A power flow supplies the true state.  Meters are generated from it with
seeded Gaussian noise, so every estimate is checked against the truth,
not just against the meters.  Then one meter at a time gets a gross
error: a bus voltage, which bad-data detection catches, and a branch
flow, which it does not.

Writes estimate.csv next to the staged inputs.
"""

import os
import random

import gridpack
from gridpack.state_estimation import Measurement
from _case import stage, stdout_to, write_csv

SEED = 14
# Meter standard deviations, per unit on 100 MVA.  Generator-bus voltages
# get the tight one: SE treats them as regulated and weights them at
# 0.001 whatever the meter says, so looser noise there reads as bad data.
SIGMA = {"VM": 0.004, "VM_gen": 0.001, "PI": 0.01, "QI": 0.01,
         "PIJ": 0.008, "QIJ": 0.008}
# (label, meter key, error in pu).  0.05 pu is 12 sigma on a voltage
# meter; 0.5 pu is 50 MW on a line carrying 73.
CASES = [
    ("clean meters", None, 0.0),
    ("VM at bus 12 reads 0.05 pu high", ("VM", 12), 0.05),
    ("PIJ 2-3 reads 0.5 pu high", ("PIJ", 2, 3), 0.5),
]
LOG = "se.log"                              # the C++ module's own output


def main():
    with gridpack.Session() as session:
        rank0 = session.rank == 0
        # One XML for both applications: Configuration merges every file
        # into one tree and reads only the first, so a second, differently
        # shaped XML would be ignored.
        study = gridpack.Case.from_data_sets("input/powerflow/input_14.xml")
        study.include("input/se/input_14.xml", "State_estimation")
        case = stage(session, study)

        # ---- 1: the true state ----------------------------------------
        with stdout_to(LOG):               # the XML turns on -ksp_view
            with gridpack.PowerFlow(session, case,
                                    suppress_output=True) as pf:
                truth = pf.solve()
                buses = truth.buses()      # collective
                branches = truth.branches()
                gens = {g["busId"] for g in truth.generators()}

        meters = make_meters(buses, branches, gens)

        # ---- 2: estimate, clean and with one bad meter ----------------
        runs = []
        for label, bad, error in CASES:
            used = [corrupt(m, error) if key(m) == bad else m for m in meters]
            runs.append((label, bad, estimate(session, case, used)))

        if rank0:
            report(session, buses, meters, runs)
            rows = []
            for b in buses:
                row = {"busId": b["busId"], "v_true": b["voltage"],
                       "angle_true": b["angle"]}
                for i, (_, _, run) in enumerate(runs):
                    row["v_case%d" % i], row["angle_case%d" % i] = (
                        run["state"][b["busId"]][1], run["state"][b["busId"]][0])
                rows.append(row)
            write_csv("estimate.csv", rows)
            print("\nwrote %s/estimate.csv (SE's own output: %s)"
                  % (os.getcwd(), LOG))


# -------------------------------------------------------------
# Meters
# -------------------------------------------------------------

def make_meters(buses, branches, gens):
    """Voltage and injection at every bus, P and Q at the sending end of
    every branch: 82 meters for 27 state variables."""
    rng = random.Random(SEED)
    meters = []

    def add(kind, value, sigma, bus=None, fbus=None, tbus=None, ckt=None):
        m = Measurement()
        m.type(kind)
        if bus is not None:
            m.busid = bus
        else:
            m.fbusid, m.tbusid = fbus, tbus
            m.ckt(ckt.strip().ljust(2)[:2])
        m.value = value + rng.gauss(0.0, sigma)
        m.deviation = sigma
        meters.append(m)

    for b in buses:
        n = b["busId"]
        add("VM", b["voltage"], SIGMA["VM_gen" if n in gens else "VM"], bus=n)
        add("PI", b["pInjection"] / 100.0, SIGMA["PI"], bus=n)
        add("QI", b["qInjection"] / 100.0, SIGMA["QI"], bus=n)
    for br in branches:
        for kind, mw in (("PIJ", br["pFrom"]), ("QIJ", br["qFrom"])):
            add(kind, mw / 100.0, SIGMA[kind], fbus=br["fromBus"],
                tbus=br["toBus"], ckt=br["circuitId"])
    return meters


def is_branch(kind):
    return kind in ("PIJ", "QIJ")


def key(m):
    if is_branch(m.type()):
        return (m.type(), m.fbusid, m.tbusid)
    return (m.type(), m.busid)


def corrupt(m, error):
    c = Measurement()
    c.type(m.type())
    if is_branch(m.type()):
        c.fbusid, c.tbusid = m.fbusid, m.tbusid
        c.ckt(m.ckt())
    else:
        c.busid = m.busid
    c.value = m.value + error
    c.deviation = m.deviation
    return c


# -------------------------------------------------------------
# Estimation
# -------------------------------------------------------------

def estimate(session, case, meters):
    """Solve SE; return the estimated state and the meters it flagged."""
    # SE prints every iteration and table with bare printf.
    with stdout_to(LOG):
        with gridpack.StateEstimation(session, case,
                                      measurements=meters) as se:
            result = se.solve()
            buses = result.buses()                # collective
            bad = result.bad_measurements()
    flagged = []
    for m in bad:
        k = ((m["type"], m["busId"]) if "busId" in m
             else (m["type"], m["fromBus"], m["toBus"]))
        if k not in flagged:
            flagged.append(k)
    return {"state": {b["busId"]: (b["angle"], b["voltage"]) for b in buses},
            "flagged": flagged, "converged": result.converged}


# -------------------------------------------------------------
# Report
# -------------------------------------------------------------

def name(k):
    return "%s %d" % k if len(k) == 2 else "%s %d-%d" % k


def report(session, buses, meters, runs):
    print("=" * 66)
    print("State estimation: IEEE 14 on %d rank(s)" % session.size)
    print("=" * 66)
    print("%d meters (seed %d) for %d state variables; the power flow is "
          "the truth" % (len(meters), SEED, 2 * len(buses) - 1))
    print("\n  %-34s %9s %10s  %s" % ("case", "max |dV|", "max |dang|",
                                     "flagged as bad"))
    for label, bad, run in runs:
        dv = max(abs(run["state"][b["busId"]][1] - b["voltage"]) for b in buses)
        da = max(abs(run["state"][b["busId"]][0] - b["angle"]) for b in buses)
        flagged = ", ".join(name(k) for k in run["flagged"]) or "none"
        print("  %-34s %7.4f pu %7.3f deg  %s" % (label, dv, da, flagged))
        if not run["converged"]:
            print("    did not converge")

    for label, bad, run in runs[1:]:
        caught = bad in run["flagged"]
        others = [k for k in run["flagged"] if k != bad]
        print("\n  %s: %s" % (label, "caught" if caught else "NOT caught"))
        if others:
            print("    also flagged, though their readings are good: %s"
                  % ", ".join(name(k) for k in others))
    # Seen on every branch-flow error tried; the ranking is r/sigma, not
    # the normalized residual, which spreads one bad meter's error.
    print("\n  No branch-flow meter is ever flagged.  The detector ranks "
          "r/sigma,\n  not the normalized residual, and here blames the "
          "injection meters\n  next to the bad flow instead.")


if __name__ == "__main__":
    main()
