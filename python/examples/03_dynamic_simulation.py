#!/usr/bin/env python
# -------------------------------------------------------------
# file: python/examples/03_dynamic_simulation.py
# -------------------------------------------------------------

"""Run the 9-bus transient case and measure the swing the fault causes.

    python 03_dynamic_simulation.py
    mpiexec -np 3 python 03_dynamic_simulation.py

Three ranks is the ceiling here: a 9-bus network cannot be partitioned
across four, and np=4 hangs in the C++ (dsf.x does the same).

This is the dsf2.x path: the XML owns the fault, GridPACK runs the whole
simulation in C++, and the per-step time series arrives as the generator
watch file.  Drive the loop from Python instead -- and pick the event at
run time -- with 04_stepped_event.py.
"""

import csv

import gridpack
from _case import stage

# From <generatorWatchFileName> in input_9b3g.xml.  <generatorWatch> picks
# the machines it covers; here that is the generator at bus 1.
WATCH_CSV = "gen_watch_9b3g.csv"

# The fault the XML declares: bus 6 - bus 7 opens at 0.03 s, clears at 0.06 s.
FAULT_BRANCH = (6, 7)


def main():
    with gridpack.Session() as session:
        case = stage(session, "input/ds/input_9b3g.xml",
                     "raw/9b3g.raw", "dyr/9b3g.dyr")

        with gridpack.DynamicSim(session, case, suppress_output=True) as ds:
            if session.rank == 0:
                print("=" * 62)
                print("Dynamic simulation: %s on %d rank(s)"
                      % (case, session.size))
                print("=" * 62)
                print("time step   : %.3f s" % ds.time_step)
                print("final time  : %.3f s" % ds.final_time)
                print("fault       : branch %d-%d\n" % FAULT_BRANCH)

            # DSFullApp::run() prints one "Time = ..." line per step with a
            # bare printf, which suppress_output cannot reach.
            result = ds.run()

            # Live queries against the final state.  These work without an
            # <observations> block, which this XML has none of; the watch
            # file below is where the time series comes from on this path.
            if session.rank == 0:
                print("\nFinal generator output (MW, MVAr):")
                for bus in (1, 2, 3):
                    print("  bus %-3d %s" % (bus, _fmt(result.get_generator_power(bus, "1"))))
                print("Final load at bus 5: %s"
                      % _fmt(result.get_bus_total_load_power(5)))

        # Outside the with-block on purpose: the CSV is GridPACK's own
        # output file and outlives the wrapper.
        if session.rank == 0:
            report_swing(WATCH_CSV)


def _fmt(pq):
    """Format a (P, Q) pair, or say so when the element was not found."""
    if pq is None:
        return "not found"
    return "%9.3f %9.3f" % (pq[0], pq[1])


def report_swing(path):
    """Summarize the watched generator's response from GridPACK's CSV."""
    with open(path) as fh:
        rows = list(csv.reader(fh))
    header = [h.strip() for h in rows[0]]
    data = [[float(x) for x in r] for r in rows[1:] if r]

    print("\nGenerator watch file: %s (%d steps, %d channels)"
          % (path, len(data), len(header) - 1))

    t = header.index("t")
    for name, unit in (("1_1_speed", "pu"), ("1_1_V", "pu"), ("1_1_Pg", "pu")):
        if name not in header:
            continue
        j = header.index(name)
        lo = min(data, key=lambda r: r[j])
        hi = max(data, key=lambda r: r[j])
        print("  %-10s min %.5f %s at t=%.2f s   max %.5f %s at t=%.2f s"
              % (name, lo[j], unit, lo[t], hi[j], unit, hi[t]))


if __name__ == "__main__":
    main()
