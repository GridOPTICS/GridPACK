#!/usr/bin/env python
# -------------------------------------------------------------
# file: python/examples/05_composed_study.py
# -------------------------------------------------------------

"""Base case, N-1 screen, dynamic check of the worst outage: one study on IEEE 14.

    python 05_composed_study.py
    mpiexec -np 2 python 05_composed_study.py

Chains three applications through one Session.  The screen's contingency
list is written from the solved base case, and the outage the screen ranks
worst is opened from Python during the dynamic run, so nothing here is
declared in an input file by hand.  Writes study.csv (the screen) and
dynamics.csv (the time series) next to the inputs.

One config file serves all three steps.  GridPACK keeps a single
process-wide Configuration tree, a later open() merges into it rather than
replacing it, and a path lookup only searches the first <Configuration>
root -- so every block any step needs has to be in the first file opened,
complete, before anything opens it.  That is also why the dynamic step
uses the stepper (demo 4's route) rather than an XML event: the event is
not known until the screen has run.
"""

import csv
import os
from collections import Counter

import gridpack
from _case import (branch_name, copy_block, copy_input, set_tag, stage,
                   write_contingency_list, write_csv)

NETWORK = "raw/IEEE14_PTIv33_rated.raw"   # all 20 branches rated
DYNAMICS = "dyr/IEEE14.dyr"
N1_LIST = "n1_branches.xml"               # written by this demo

# Dynamic check: open the worst branch at TRIP_TIME and watch every machine
# until T_END.  The stock rtpr input runs 1 s, too short to see whether the
# swing settles.
TRIP_TIME = 1.0
T_END = 5.0

# Finalist order when nothing lifts a bus out of the band.  Cases the
# power flow cannot solve are reported but not simulated: there is no
# post-contingency operating point to start from.
_SEVERITY = {"BUS+BRANCH VIOLATION": 0, "BUS VIOLATION": 1,
             "BRANCH VIOLATION": 2}


def main():
    with gridpack.Session() as session:
        rank0 = session.rank == 0
        case = stage(session, "input/ca/input_14.xml", network=NETWORK)
        dyr = copy_input(session, DYNAMICS, os.path.basename(DYNAMICS))

        # ---- 0: one config for the whole study --------------------------
        # The CA input supplies Contingency_analysis and Powerflow; the
        # rtpr input supplies a Dynamic_simulation block already paired
        # with IEEE14.dyr.  The horizon is extended and every machine
        # observed.  Its XML fault is kept but pushed past the horizon, not
        # deleted: the stepper's initializer takes the XML event list when
        # Python passes none and reads its first entry unconditionally
        # (hadrec_app_module.cpp, initializeDynSimu), so an empty list is a
        # segfault.  The event this study wants is applied from Python.
        if rank0:
            set_tag(case, "contingencyList", N1_LIST,
                    parent="Contingency_analysis")
            copy_block("input/rtpr/input_14.xml", case, "Dynamic_simulation")
            set_tag(case, "simulationTime", T_END)
            set_tag(case, "beginFault", 2 * T_END)
            set_tag(case, "endFault", 2 * T_END)
            set_tag(case, "observations",
                    observations(read_generators(dyr),
                                 read_buses(os.path.basename(NETWORK))),
                    parent="Dynamic_simulation")
        session.barrier()

        # ---- 1: base case --------------------------------------------
        with gridpack.PowerFlow(session, case, suppress_output=True) as pf:
            base = pf.solve()
            # Collective: every rank calls them.
            buses = base.buses()
            branches = base.branches()
        if rank0:
            report_base(case, session, base, buses, branches)
            # The XML asks for FullBranchN1, which pygridpack does not
            # expand; the list is read lazily, so writing it now is in time.
            write_contingency_list(N1_LIST, branches)
        session.barrier()

        # ---- 2: N-1 screen, 3: pick the worst --------------------------
        with gridpack.ContingencyAnalysis(session, case,
                                          print_calc_files=False,
                                          suppress_output=True) as ca:
            ca.run()
            results = ca.gather()

            # The screen says whether a limit was crossed, not by how much:
            # re-solve the finalists to rank them.
            finalists = sorted((r for r in results if r.status in _SEVERITY),
                               key=lambda r: _SEVERITY[r.status])
            measured = [measure(ca, r) for r in finalists]
            worst = pick_worst(measured)
            vmin, vmax = ca.min_voltage, ca.max_voltage

        if rank0:
            report_screen(results, measured, worst, vmin, vmax)
            write_csv("study.csv", study_rows(results, measured))

        if worst is None:
            if rank0:
                print("\nNo solvable outage violates a limit; nothing to "
                      "verify dynamically.\nwrote %s/study.csv" % os.getcwd())
            return

        # ---- 4: open the worst outage in a dynamic run -----------------
        with gridpack.DynamicSimStepper(session, case,
                                        suppress_output=True) as st:
            if rank0:
                print("\n" + "=" * 66)
                print("Dynamic check: open %s at t = %.1f s, run to %.1f s"
                      % (worst["name"], TRIP_TIME, T_END))
                print("=" * 66)
            tripped = False
            while not st.done:
                st.step()
                if not tripped and st.current_time >= TRIP_TIME:
                    st.apply_line_trip(worst["from"], worst["to"],
                                       ckt=worst["ckt"])
                    tripped = True
            # Rank-local lists that every rank holds alike, so one writer.
            if rank0:
                report_recovery(st.result.to_records(),
                                st.result.channel_names(), worst)
                st.result.to_csv("dynamics.csv")
                print("\nwrote %s/{study,dynamics}.csv" % os.getcwd())


# -------------------------------------------------------------
# step 0 helpers
# -------------------------------------------------------------

def read_generators(dyr_path):
    """(bus, id) for each machine in a .dyr: `bus, 'MODEL', 'id', ...`."""
    gens = []
    for line in open(dyr_path):
        parts = line.split(",")
        if len(parts) > 2:
            gens.append((int(parts[0]), parts[2].strip().strip("'").strip()))
    return gens


def read_buses(raw_path):
    """Bus numbers from a PSS/E RAW file: the first block after the header."""
    buses = []
    with open(raw_path) as fh:
        for line in list(fh)[3:]:
            if line.startswith("0 "):
                break
            buses.append(int(line.split(",")[0]))
    return buses


def observations(gens, buses):
    """An <observations> body: every machine and every bus."""
    obs = ["      <observation><type>generator</type><busID>%d</busID>"
           "<generatorID>%s</generatorID></observation>" % g for g in gens]
    obs += ["      <observation><type>bus</type><busID>%d</busID>"
            "</observation>" % b for b in buses]
    return "\n" + "\n".join(obs) + "\n    "


# -------------------------------------------------------------
# step 3 helpers
# -------------------------------------------------------------

def measure(ca, r):
    """Re-solve one screened contingency and read its worst numbers.

    Same sequence ContingencyAnalysis uses per case, on the same private
    network, so every rank computes the same answer without a broadcast.
    """
    pf = ca.powerflow
    ctg = next(c for c in ca.contingencies if c.name == r.name)
    pyb = ctg.to_pybind()
    pf.reset_voltages()
    pf.set_contingency(pyb)
    try:
        res = pf.solve(strict=False)
        if res.converged and ca.check_qlim and not pf.check_qlim_violations():
            res = pf.solve(strict=False)
        v = res.violations(min_voltage=ca.min_voltage,
                           max_voltage=ca.max_voltage)
        loading = max((b["loadingPercent"] for b in res.branches()
                       if b["rateA"] > 0), default=0.0)
        low = min(res.buses(), key=lambda b: b["voltage"])
    finally:
        pf.unset_contingency(pyb)
        if ca.check_qlim:
            pf.clear_qlim_violations()
    return {"name": r.name, "status": r.status,
            "from": ctg.from_buses[0], "to": ctg.to_buses[0],
            "ckt": ctg.circuit_ids[0].strip(),
            "bus_violations": len(v["voltage"]),
            "low_bus": low["busId"], "low_name": low["name"],
            "low_voltage": low["voltage"],
            "max_loading": loading}


def pick_worst(measured):
    """The outage that pulls a bus out of the band, else the largest overload.

    Ties on bus count break on the lowest voltage.  This criterion is the
    demo's, not GridPACK's.
    """
    if not measured:
        return None
    bus = [m for m in measured if m["bus_violations"]]
    if bus:
        return min(bus, key=lambda m: (-m["bus_violations"], m["low_voltage"]))
    return max(measured, key=lambda m: m["max_loading"])


def study_rows(results, measured):
    by_name = {m["name"]: m for m in measured}
    rows = []
    for r in results:
        m = by_name.get(r.name, {})
        rows.append({"contingency": r.name.strip(), "status": r.status,
                     "bus_violations": m.get("bus_violations", ""),
                     "low_bus": m.get("low_bus", ""),
                     "low_voltage": m.get("low_voltage", ""),
                     "max_loading_pct": m.get("max_loading", "")})
    return rows


# -------------------------------------------------------------
# reports
# -------------------------------------------------------------

def report_base(case, session, base, buses, branches):
    print("=" * 66)
    print("Composed study: %s on %d rank(s)" % (case, session.size))
    print("=" * 66)
    print("base case   : converged in %s iterations, %d buses, %d branches"
          % (base.iterations, len(buses), len(branches)))

    print("\nLowest voltages:")
    for b in sorted(buses, key=lambda b: b["voltage"])[:3]:
        print("  bus %-3d %-8s %.4f pu" % (b["busId"], b["name"], b["voltage"]))

    print("Most loaded branches:")
    hot = sorted((b for b in branches if b["rateA"] > 0),
                 key=lambda b: -b["loadingPercent"])[:3]
    for b in hot:
        print("  %-10s %5.1f%% of %g MVA"
              % (branch_name(b), b["loadingPercent"], b["rateA"]))


def report_screen(results, measured, worst, vmin, vmax):
    print("\n" + "=" * 66)
    print("N-1 screen: %d single-branch outages, band [%.2f, %.2f] pu"
          % (len(results), vmin, vmax))
    print("=" * 66)
    for status, n in sorted(Counter(r.status for r in results).items(),
                            key=lambda kv: -kv[1]):
        print("  %-24s %3d" % (status, n))

    unsolved = [r for r in results if r.status not in _SEVERITY
                and r.status != "OK"]
    if unsolved:
        print("\nNot simulated (no post-contingency operating point):")
        for r in unsolved:
            print("  %-10s %s" % (r.name.strip(), r.status))

    if measured:
        print("\nFinalists, re-solved:")
        print("  %-10s %-22s %5s  %-19s %8s"
              % ("outage", "status", "buses", "lowest bus", "loading"))
        for m in measured:
            print("  %-10s %-22s %5d  %-3d %-8s %.4f %7.1f%%"
                  % (m["name"], m["status"], m["bus_violations"],
                     m["low_bus"], m["low_name"], m["low_voltage"],
                     m["max_loading"]))
    if worst:
        why = ("%d bus(es) out of band" % worst["bus_violations"]
               if worst["bus_violations"]
               else "largest overload, %.1f%%" % worst["max_loading"])
        print("\nWorst: %s (%s)" % (worst["name"], why))


def report_recovery(records, channels, worst):
    """What the trip did, from the observation series.

    Voltages at generator buses come back constant at their setpoint on
    this path, so the trip is shown on machine power and the voltage
    survey is worth reading only for load buses.
    """
    import math
    before = [r for r in records if r["time"] < TRIP_TIME][-1]
    # The trip is applied after the step that reaches TRIP_TIME, so its
    # first visible effect is one step later.
    after = [r for r in records if r["time"] > TRIP_TIME]
    t_end = after[-1]["time"]
    machines = [c[4:-7] for c in channels if c.endswith("_rspeed")]
    print("\nObservations: %d steps, %d channels, %d after the trip"
          % (len(records), len(channels) - 1, len(after)))

    # The trip's footprint: the machine whose output jumps most across
    # the event.  No jump means nothing opened.
    dp, m = max((abs(after[0]["gen_%s_P" % m] - before["gen_%s_P" % m]), m)
                for m in machines)
    print("  trip footprint: gen %s P %.3f -> %.3f pu at t = %.2f s"
          % (m, before["gen_%s_P" % m], after[0]["gen_%s_P" % m],
             after[0]["time"]))

    print("\n  Lowest voltage after the trip:")
    lows = sorted(((min(r[c] for r in after), c[4:-5])
                   for c in channels if c.endswith("_vmag")))[:3]
    for v, bus in lows:
        print("  bus %-4s %.4f pu" % (bus, v))

    # Swing amplitude in the first second after the trip against the
    # last: near-equal is an undamped swing.  Angles are read relative to
    # the first machine: without governors the whole system's speed
    # drifts, which moves every absolute angle alike and says nothing
    # about synchronism.  A relative angle that has moved more than pi
    # is a machine that slipped poles.
    early = [r for r in after if r["time"] < TRIP_TIME + 1.0]
    late = [r for r in after if r["time"] >= t_end - 1.0]
    ref = "gen_%s_rangle" % machines[0]
    print("\n  %-8s %9s %9s %11s %11s %14s"
          % ("machine", "min spd", "max spd", "swing 1st s", "swing last",
             "angle vs " + machines[0]))
    slipped = []
    for m in machines:
        spd, ang = "gen_%s_rspeed" % m, "gen_%s_rangle" % m
        s = [r[spd] for r in after]
        e = [r[spd] for r in early]
        l = [r[spd] for r in late]
        rel = ((after[-1][ang] - after[-1][ref])
               - (before[ang] - before[ref]))
        if abs(rel) > math.pi:
            slipped.append(m)
        print("  %-8s %9.5f %9.5f %11.5f %11.5f %10.3f rad"
              % (m, min(s), max(s), max(e) - min(e), max(l) - min(l), rel))

    mean_spd = sum(after[-1]["gen_%s_rspeed" % m] for m in machines) / len(machines)
    if slipped:
        print("\n  Verdict: %s lost synchronism by %.1f s; the outage is not "
              "dynamically secure." % (", ".join("gen " + m for m in slipped),
                                       t_end))
    else:
        print("\n  Verdict: all machines in step at %.1f s; system speed "
              "%.4f pu (%.2f Hz), still drifting with no governors."
              % (t_end, mean_spd, 60.0 * mean_spd))


if __name__ == "__main__":
    main()
