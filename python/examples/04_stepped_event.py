#!/usr/bin/env python
# -------------------------------------------------------------
# file: python/examples/04_stepped_event.py
# -------------------------------------------------------------

"""Trip a generator from Python mid-simulation and measure the frequency dip.

    python 04_stepped_event.py

Writes stepped_event.csv next to the staged inputs.

The step loop lives in Python, so the event is decided at run time rather
than declared in the XML -- this is the shape an RL agent or a co-simulation
driver needs.  For a whole run handed to C++ in one call, see
03_dynamic_simulation.py.
"""

import os

import gridpack
from _case import stage

# Kundur's two-area case: 4 machines, 2 areas joined by a weak tie.  Its
# <observations> block is what makes the per-step channels available;
# without one the stepper collects nothing.
CASE = "input/ds/input_kundur_two_area_wsieg1.xml"

TRIP_BUS, TRIP_GEN = 2, "1"      # 700 MW machine in area 1
TRIP_TIME = 1.0                  # s, after the flat start settles
T_END = 8.0                      # s; the XML declares 20 s and a fault at 10 s


def main():
    with gridpack.Session() as session:
        case = stage(session, CASE, "raw/kundur-twoarea_v33.raw",
                     "dyr/kundur-twoarea_wsieg1.dyr")

        with gridpack.DynamicSimStepper(session, case,
                                        suppress_output=True) as st:
            if session.rank == 0:
                print("=" * 62)
                print("Stepped run: %s on %d rank(s)" % (case, session.size))
                print("=" * 62)
                print("time step   : %.4f s" % st.time_step)
                print("stopping at : %.1f s" % T_END)
                print("event       : trip generator %s at bus %d, t = %.1f s\n"
                      % (TRIP_GEN, TRIP_BUS, TRIP_TIME))

            tripped = False
            # Stop on our own clock, not on st.done: the XML runs to 20 s and
            # faults at 10 s, and the injected trip is meant to be the only
            # disturbance in the window.
            while not st.done and st.current_time < T_END:
                st.step()
                if not tripped and st.current_time >= TRIP_TIME:
                    st.apply_generator_trip(bus_number=TRIP_BUS,
                                            gen_id=TRIP_GEN)
                    tripped = True

            # Per-rank local lists, not a gather: every rank holds the same
            # series, so one writer is enough and guarding is safe here
            # (unlike the power-flow results API, which gathers inside).
            if session.rank == 0:
                report(st)
                st.result.to_csv("stepped_event.csv")
                print("\nwrote %s/stepped_event.csv" % os.getcwd())


def report(st):
    records = st.result.to_records()
    channels = st.result.channel_names()
    print("steps recorded: %d, last t = %.3f s" % (len(records),
                                                   records[-1]["time"]))
    print("channels      : %s\n" % ", ".join(channels[1:]))

    def series(name):
        return [(r["time"], r[name]) for r in records]

    # Frequency: how far it fell, and whether the governors caught it before
    # the window closed.
    for name in [c for c in channels if c.endswith("_rspeed")]:
        s = series(name)
        t_low, low = min(s, key=lambda x: x[1])
        trend = "recovering" if s[-1][1] > low else "still falling"
        print("%-16s nadir %.5f pu at t=%.2f s, %s (%.5f pu at %.1f s)"
              % (name, low, t_low, trend, s[-1][1], s[-1][0]))

    print()
    for name in [c for c in channels if c.endswith("_vmag")]:
        s = series(name)
        t_low, low = min(s, key=lambda x: x[1])
        print("%-16s dips to %.4f pu at t=%.2f s, %.4f pu at the end"
              % (name, low, t_low, s[-1][1]))

    print()
    for name in [c for c in channels if c.endswith("_P")]:
        s = series(name)
        print("%-16s %.3f -> %.3f pu (peak %.3f)"
              % (name, s[0][1], s[-1][1], max(v for _, v in s)))


if __name__ == "__main__":
    main()
