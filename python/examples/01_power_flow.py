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

VMIN, VMAX = 0.95, 1.05          # tighter than the 0.9/1.1 default, to show a hit


def main():
    # One Session per process, and it owns shutdown ordering: build it first
    # and let the with-block close it, or MPI_Finalize runs before GridPACK
    # has let go of the network.
    with gridpack.Session() as session:
        case = stage(session, "input/powerflow/input_14.xml", "raw/IEEE14.raw")

        pf = gridpack.PowerFlow(session, case, suppress_output=True)

        # strict=False returns the unconverged result instead of raising, so
        # the report below can say what went wrong.  Drop it to get a
        # PowerFlowDiverged exception on every rank instead.
        result = pf.solve(strict=False)

        # Collective: these gather across ranks, so every rank has to reach
        # them.  Guarding them with `if rank == 0` deadlocks.
        buses = result.buses()
        violations = result.violations(min_voltage=VMIN, max_voltage=VMAX)

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
        for v in violations["voltage"]:
            print("  bus %-6d %.4f pu (%s)" % (v["busId"], v["voltage"], v["kind"]))
    else:
        print("  none")

    print("\nBranches loaded over %.0f%%:" % band["overload_threshold"])
    if violations["overload"]:
        for o in violations["overload"]:
            print("  %d->%-4d %-3s %6.1f%% of %.0f MVA"
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
