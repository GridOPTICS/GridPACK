#!/usr/bin/env python
# -------------------------------------------------------------
# file: python/examples/01_power_flow.py
# -------------------------------------------------------------

"""Solve a power flow on IEEE 14 and report what is outside limits.

    python 01_power_flow.py
    mpiexec -np 4 python 01_power_flow.py

Writes buses.csv next to the staged inputs.
"""

import os

import gridpack
from _case import stage, write_csv

# Tighter than the 0.9/1.1 and 100% defaults, so the base case has hits:
# five branches sit between 76% and 84% loaded with nothing switched out.
VMIN, VMAX = 0.95, 1.05
OVERLOAD = 75.0


def main():
    # One Session per process, and it owns shutdown ordering: build it first
    # and let the with-block close it, or MPI_Finalize runs before GridPACK
    # has let go of the network.
    with gridpack.Session() as session:
        # IEEE14.raw, which this XML names, carries no branch ratings, so
        # the overload check would have nothing to test.  This v33 variant
        # rates all 20; the version is auto-detected from the file header.
        case = stage(session, "input/powerflow/input_14.xml",
                     network="raw/IEEE14_PTIv33_rated.raw")

        pf = gridpack.PowerFlow(session, case, suppress_output=True)

        # strict=False returns the unconverged result instead of raising, so
        # the report below can say what went wrong.  Drop it to get a
        # PowerFlowDiverged exception on every rank instead.
        result = pf.solve(strict=False)

        # Collective: these gather across ranks, so every rank has to reach
        # them.  Guarding them with `if rank == 0` deadlocks.
        buses = result.buses()
        violations = result.violations(min_voltage=VMIN, max_voltage=VMAX,
                                       overload_threshold=OVERLOAD)

        if session.rank == 0:
            report(case, session, result, buses, violations)
            write_csv("buses.csv", buses)
            print("\nwrote %s/buses.csv" % os.getcwd())


def report(case, session, result, buses, violations):
    print("=" * 62)
    print("Power flow: %s on %d rank(s)" % (case, session.size))
    print("=" * 62)
    print("converged   : %s" % result.converged)
    print("iterations  : %s" % result.iterations)
    print("buses       : %d" % len(buses))
    print("branches    : %d" % violations["n_branches"])

    if not result.converged:
        print("\nworst mismatch: %s" % (result.mismatch,))
        return

    band = violations["limits"]
    print("\nVoltage outside [%.2f, %.2f] pu:" % (band["min_voltage"],
                                                  band["max_voltage"]))
    if violations["voltage"]:
        # The name is the RAW file's, so it is whatever the model author
        # wrote: generic here, real substation names on IEEE 118 (demo 2).
        for v in violations["voltage"]:
            print("  bus %-4d %-12s %.4f pu (%s)"
                  % (v["busId"], v["name"], v["voltage"], v["kind"]))
    else:
        print("  none")

    print("\nBranches loaded over %.0f%%:" % band["overload_threshold"])
    if violations["overload"]:
        for o in violations["overload"]:
            print("  %d->%-4d %-3s %6.1f%% of %g MVA"
                  % (o["fromBus"], o["toBus"], o["circuitId"],
                     o["loadingPercent"], o["rateA"]))
    else:
        print("  none")

    # Without this an unrated network looks like a network with no overloads.
    if violations["unrated_branches"]:
        print("  (%d of %d branches carry no rating and were not checked)"
              % (violations["unrated_branches"], violations["n_branches"]))


if __name__ == "__main__":
    main()
