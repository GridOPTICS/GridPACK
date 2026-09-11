#!/usr/bin/env python
# -------------------------------------------------------------
# file: python/examples/02_contingency_screening.py
# -------------------------------------------------------------

"""Screen 179 N-1 line outages on IEEE 118 and rank what they break.

    python 02_contingency_screening.py
    mpiexec -np 4 python 02_contingency_screening.py

Contingencies are handed out across ranks, so more ranks finishes the same
screen sooner.  Writes screening.csv next to the staged inputs.
"""

import os
from collections import Counter

import gridpack
from _case import stage, write_csv


def main():
    with gridpack.Session() as session:
        # The XML names contingencies_118.xml and IEEE118.raw as bare
        # filenames, so all three have to land in one directory.
        case = stage(session,
                     "input/ca/input_118.xml",
                     "contingencies/contingencies_118.xml",
                     "raw/IEEE118.raw")

        # print_calc_files=False: 179 per-contingency .out files are not
        # what this demo is showing.
        with gridpack.ContingencyAnalysis(session, case,
                                          print_calc_files=False,
                                          suppress_output=True) as ca:
            if not ca.contingencies:
                # pygridpack needs an explicit <contingencyList>; it cannot
                # expand FullBranchN1/FullGeneratorN1 yet.
                raise SystemExit(
                    "no contingencies in %s -- set <contingencyList>" % case)

            if session.rank == 0:
                print("=" * 66)
                print("N-1 screen: %d contingencies on %d rank(s)"
                      % (len(ca.contingencies), session.size))
                print("band [%.2f, %.2f] pu"
                      % (ca.min_voltage, ca.max_voltage))
                print("=" * 66)

            base = ca.solve_base_case()
            if not base.converged:
                raise SystemExit("base case did not converge; nothing to "
                                 "measure contingencies against")

            # Collective, so every rank calls it.  Worth printing: only
            # voltage violations are exempted below, never overloads, so a
            # base case that is already overloaded flags every contingency.
            pre = base.violations(min_voltage=ca.min_voltage,
                                  max_voltage=ca.max_voltage)
            if session.rank == 0:
                print("base case: %d bus / %d branch violations already "
                      "present" % (len(pre["voltage"]), len(pre["overload"])))
                print("  (bus violations here are exempted from the screen; "
                      "overloads are not)\n")

            ca.run()

            # Every rank evaluated a disjoint subset; gather is collective.
            results = ca.gather()

        if session.rank == 0:
            summarise(results)
            write_csv("screening.csv", [
                {"name": r.name.strip(), "status": r.status,
                 "converged": r.converged, "voltage_ok": r.voltage_ok,
                 "overload_ok": r.overload_ok, "islanded": r.islanded,
                 "slack_overload": r.slack_overload}
                for r in results])
            print("\nwrote %s/screening.csv" % os.getcwd())


def summarise(results):
    print("Outcome                    Cases")
    print("-" * 40)
    for status, n in sorted(Counter(r.status for r in results).items(),
                            key=lambda kv: -kv[1]):
        print("  %-24s %4d" % (status, n))

    # The cases worth a human's time: not solvable at all, or a bus outside
    # the band that the base case was not already outside.
    unsolved = [r for r in results if not r.converged]
    bus = [r for r in results if r.converged and not r.voltage_ok]

    print("\n%d of %d could not be solved:" % (len(unsolved), len(results)))
    for r in unsolved:
        print("  %-10s %s" % (r.name.strip(), r.report_lines[0]))

    print("\n%d introduced a new bus violation:" % len(bus))
    for r in bus:
        print("  %-10s %s" % (r.name.strip(), r.status))


if __name__ == "__main__":
    main()
