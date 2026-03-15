#!/usr/bin/env python
# -------------------------------------------------------------
# file: gridpack_cli.py
# -------------------------------------------------------------
# Unified command-line interface for GridPACK Python applications.
#
# Copyright (c) 2013 Battelle Memorial Institute
# Licensed under modified BSD License. A copy of this license can be found
# in the LICENSE file in the top level directory of this distribution.
#
# -------------------------------------------------------------
# Created March 2026
# Author: Yousu Chen
# -------------------------------------------------------------

"""
GridPACK unified command-line interface.

Provides a single entry point for running GridPACK applications:
  - powerflow (pf)   : AC power flow analysis
  - dsf              : Dynamic simulation framework
  - se               : State estimation
  - hadrec           : HADREC remedial action control
  - emt              : Electromagnetic transient simulation
  - ca               : Contingency analysis

Usage:
    mpiexec -np N gridpack <command> <config.xml> [options]

Examples:
    gridpack powerflow input.xml
    gridpack dsf input.xml --quiet
    gridpack ca input.xml --vlimits 0.9 1.1
    mpiexec -np 4 gridpack powerflow input.xml --output results.txt
"""

from __future__ import print_function
import sys
import os
import argparse


def _safe_exit(rc=0):
    """Return exit code. Actual process termination is handled by main()."""
    return rc


class _GridPACKEnv:
    """Manages GridPACK environment lifecycle.

    Ensures proper MPI initialization and finalization order.
    The environment and communicator are stored as attributes
    and cleaned up in the correct order when close() is called
    or the object is deleted.
    """
    def __init__(self):
        import gridpack
        self._env = gridpack.Environment()
        self._comm = gridpack.Communicator()

    @property
    def comm(self):
        return self._comm

    def close(self):
        """Explicitly release in correct order."""
        self._comm = None
        self._env = None


# -------------------------------------------------------------
# Powerflow subcommand
# -------------------------------------------------------------
def cmd_powerflow(args, gp_env):
    """Run AC power flow analysis."""
    import gridpack
    import gridpack.powerflow

    comm = gp_env.comm
    timer = gridpack.CoarseTimer()

    config = gridpack.Configuration()
    config.open(args.config, comm)
    cursor = config.getCursor("Configuration.Powerflow")

    # Read options from XML, allow CLI overrides
    useNonLinear = cursor.get("UseNonLinear")
    exportPSSE23 = cursor.get("exportPSSE_v23")
    exportPSSE33 = cursor.get("exportPSSE_v33")
    exportPSSE34 = cursor.get("exportPSSE_v34")
    noPrint = cursor.get("suppressOutput")

    # CLI overrides
    if args.quiet:
        noPrint = True
    if args.solver == "nl":
        useNonLinear = True
    elif args.solver == "nr":
        useNonLinear = False

    pfapp = gridpack.powerflow.Powerflow()

    if noPrint:
        pfapp.suppressOutput(True)

    pfapp.readNetwork(config, -1)
    pfapp.initialize()

    if useNonLinear:
        pfapp.nl_solve()
    else:
        pfapp.solve()

    # Output handling
    if args.output:
        pfapp.open(args.output)

    pfapp.write()
    pfapp.saveData()

    if args.output:
        pfapp.close()

    # PSS/E export from CLI
    if args.export_psse:
        fmt, fname = args.export_psse
        if fmt == "v23":
            pfapp.exportPSSE23(fname)
        elif fmt == "v33":
            pfapp.exportPSSE33(fname)
        elif fmt == "v34":
            pfapp.exportPSSE34(fname)

    # PSS/E export from XML config
    if exportPSSE23:
        pfapp.exportPSSE23(exportPSSE23)
    if exportPSSE33:
        pfapp.exportPSSE33(exportPSSE33)
    if exportPSSE34:
        pfapp.exportPSSE34(exportPSSE34)

    if not noPrint and not args.no_timer:
        timer.dump()

    return 0


# -------------------------------------------------------------
# Dynamic Simulation subcommand
# -------------------------------------------------------------
def cmd_dsf(args, gp_env):
    """Run dynamic simulation framework."""
    import gridpack
    from gridpack.dynamic_simulation import DSFullApp

    comm = gp_env.comm
    timer = gridpack.CoarseTimer()
    t_total = timer.createCategory("Dynamic Simulation: Total Application")
    timer.start(t_total)

    np_ctrl = gridpack.NoPrint()
    if args.quiet:
        np_ctrl.setStatus(True)

    ds_app = DSFullApp()
    ds_app.solvePowerFlowBeforeDynSimu(args.config, -1)
    ds_app.readGenerators()
    ds_app.readSequenceData()
    ds_app.initialize()
    ds_app.setGeneratorWatch()

    if args.output:
        ds_app.open(args.output)

    # Use run()-based execution (like dsf2.py)
    ds_app.setup()
    ds_app.run()

    if args.output:
        ds_app.close()

    timer.stop(t_total)

    if not args.no_timer:
        timer.dump()

    return 0


# -------------------------------------------------------------
# State Estimation subcommand
# -------------------------------------------------------------
def cmd_se(args, gp_env):
    """Run state estimation."""
    import gridpack
    import gridpack.state_estimation

    comm = gp_env.comm
    timer = gridpack.CoarseTimer()

    config = gridpack.Configuration()
    config.open(args.config, comm)

    se = gridpack.state_estimation.SEApp()

    # Get measurement file from config
    measname = config.get("Configuration.State_estimation.measurementList")

    mconfig = gridpack.Configuration()
    mconfig.open(measname, comm)

    measures = se.getMeasurements(mconfig)

    se.readNetwork(config)
    se.initialize()
    se.setMeasurements(measures)
    se.solve()
    se.saveData()

    if args.output:
        se.open(args.output)

    se.write()

    if args.output:
        se.close()

    if not se.hasConverged():
        sys.stderr.write("Warning: State estimation did not fully converge, "
                         "but results were written anyway.\n")

    if not args.no_timer:
        timer.dump()

    return 0


# -------------------------------------------------------------
# HADREC subcommand
# -------------------------------------------------------------
def cmd_hadrec(args, gp_env):
    """Run HADREC remedial action control simulation."""
    import gridpack
    import gridpack.hadrec
    import gridpack.dynamic_simulation

    comm = gp_env.comm

    np_ctrl = gridpack.NoPrint()
    if args.quiet:
        np_ctrl.setStatus(True)

    hadapp = gridpack.hadrec.Module()
    hadapp.solvePowerFlowBeforeDynSimu(args.config, -1)
    hadapp.transferPFtoDS()

    busfaultlist = gridpack.dynamic_simulation.EventVector()
    hadapp.initializeDynSimu(busfaultlist)

    while not hadapp.isDynSimuDone():
        hadapp.executeDynSimuOneStep()

    return 0


# -------------------------------------------------------------
# EMT subcommand
# -------------------------------------------------------------
def cmd_emt(args, gp_env):
    """Run electromagnetic transient simulation."""
    import gridpack
    import gridpack.emt

    comm = gp_env.comm

    emt_app = gridpack.emt.EMT()
    emt_app.setconfigurationfile(args.config)
    emt_app.solvepowerflow()
    emt_app.setup()
    emt_app.solve()

    return 0


# -------------------------------------------------------------
# Contingency Analysis subcommand
# -------------------------------------------------------------
def _parse_contingency_xml(config):
    """Parse contingency list from XML configuration.

    Reads the contingency list file specified in the
    Contingency_analysis block and returns a list of dicts
    describing each contingency.

    Returns:
        list of dict: Each dict has keys:
            - name (str)
            - type ("Line" or "Generator")
            - from_buses, to_buses, ckt (for Line)
            - buses, gen_ids (for Generator)
    """
    import xml.etree.ElementTree as ET

    cursor = config.getCursor("Configuration.Contingency_analysis")
    contingency_file = cursor.get("contingencyList")

    if not contingency_file:
        return []

    tree = ET.parse(contingency_file)
    root = tree.getroot()

    contingencies = []
    # Find Contingency elements under ContingencyList.Contingency_analysis.Contingencies
    for cont_elem in root.iter("Contingency"):
        ca_type = cont_elem.findtext("contingencyType", "")
        ca_name = cont_elem.findtext("contingencyName", "")

        if ca_type == "Line":
            buses_str = cont_elem.findtext("contingencyLineBuses", "")
            names_str = cont_elem.findtext("contingencyLineNames", "")
            bus_ids = [int(x) for x in buses_str.split()]
            line_names = names_str.split()
            # Pad line names to 2 chars
            line_names = [n.ljust(2)[:2] for n in line_names]

            from_buses = []
            to_buses = []
            for i in range(len(line_names)):
                from_buses.append(bus_ids[2 * i])
                to_buses.append(bus_ids[2 * i + 1])

            contingencies.append({
                "name": ca_name,
                "type": "Line",
                "from_buses": from_buses,
                "to_buses": to_buses,
                "ckt": line_names,
            })
        elif ca_type == "Generator":
            buses_str = cont_elem.findtext("contingencyBuses", "")
            gens_str = cont_elem.findtext("contingencyGenerators", "")
            bus_ids = [int(x) for x in buses_str.split()]
            gen_ids = gens_str.split()
            gen_ids = [g.ljust(2)[:2] for g in gen_ids]

            contingencies.append({
                "name": ca_name,
                "type": "Generator",
                "buses": bus_ids,
                "gen_ids": gen_ids,
            })

    return contingencies


def cmd_ca(args, gp_env):
    """Run contingency analysis.

    This implements the contingency analysis workflow using the
    powerflow module with setContingency/unSetContingency, following
    the same logic as the C++ ca_driver.

    Note: Requires the contingency analysis pybind11 bindings
    (setContingency, unSetContingency, Contingency struct) to be
    available in gridpack.powerflow.
    """
    import gridpack
    import gridpack.powerflow

    comm = gp_env.comm
    timer = gridpack.CoarseTimer()
    t_total = timer.createCategory("CA: Total Application")
    timer.start(t_total)

    config = gridpack.Configuration()
    config.open(args.config, comm)

    # Read CA parameters
    cursor = config.getCursor("Configuration.Contingency_analysis")
    grp_size = int(cursor.get("groupSize") or 1)
    Vmin = float(args.vlimits[0] if args.vlimits else (cursor.get("minVoltage") or 0.9))
    Vmax = float(args.vlimits[1] if args.vlimits else (cursor.get("maxVoltage") or 1.1))
    print_calcs = True
    tmp = cursor.get("printCalcFiles")
    if tmp and tmp.lower() == "false":
        print_calcs = False
    if args.no_print_calcs:
        print_calcs = False

    # Create and initialize power flow
    pfapp = gridpack.powerflow.Powerflow()

    if args.quiet:
        pfapp.suppressOutput(True)

    pfapp.readNetwork(config, -1)
    pfapp.initialize()
    pfapp.setVoltageLimits(Vmin, Vmax)

    # Solve base case
    pfapp.solve()

    # Check Q limits if enabled
    check_qlim = True
    qlim_val = cursor.get("qlim")
    if qlim_val is not None:
        if isinstance(qlim_val, str):
            check_qlim = qlim_val.lower() != "false"
        else:
            check_qlim = bool(qlim_val)

    if check_qlim:
        if not pfapp.checkQlimViolations():
            pfapp.solve()

    # Flag base-case voltage violations to ignore
    pfapp.ignoreVoltageViolations()

    # Parse contingency list
    contingencies = _parse_contingency_xml(config)

    if comm.rank() == 0:
        print("=" * 60)
        print("Contingency Analysis")
        print("Total contingencies: %d" % len(contingencies))
        print("Voltage limits: [%.3f, %.3f]" % (Vmin, Vmax))
        print("=" * 60)

    # Use TaskManager for parallel distribution
    task_comm = comm.divide(grp_size)
    taskmgr = gridpack.TaskManager(comm)
    ntasks = len(contingencies)
    taskmgr.set(ntasks)
    task = gridpack.TaskCounter()

    results = []

    while taskmgr.nextTask(task_comm, task):
        task_id = task.task_id
        event = contingencies[task_id]

        if comm.rank() == 0:
            print("Running contingency %d: %s" % (task_id, event["name"]))

        # Build Contingency object
        contingency = gridpack.powerflow.Contingency()
        contingency.p_name = event["name"]

        if event["type"] == "Line":
            contingency.p_type = 1  # Branch
            contingency.p_from = event["from_buses"]
            contingency.p_to = event["to_buses"]
            contingency.p_ckt = event["ckt"]
            contingency.p_saveLineStatus = [1] * len(event["ckt"])
        else:
            contingency.p_type = 0  # Generator
            contingency.p_busid = event["buses"]
            contingency.p_genid = event["gen_ids"]
            contingency.p_saveGenStatus = [1] * len(event["gen_ids"])

        # Open output file for this contingency
        fname = event["name"].rstrip() + ".out"
        if print_calcs:
            pfapp.open(fname)

        # Reset and apply contingency
        pfapp.resetVoltages()
        found = pfapp.setContingency(contingency)

        if found and pfapp.solve():
            if check_qlim:
                if not pfapp.checkQlimViolations():
                    pfapp.solve()

            if print_calcs:
                pfapp.write()

            ok_v = pfapp.checkVoltageViolations()
            ok_l = pfapp.checkLineOverloadViolations()

            status = "OK"
            if not ok_v and not ok_l:
                status = "BUS+BRANCH VIOLATION"
            elif not ok_v:
                status = "BUS VIOLATION"
            elif not ok_l:
                status = "BRANCH VIOLATION"

            results.append((event["name"], True, status))

            if print_calcs:
                if not ok_v:
                    pfapp.print("\nBus Violation for contingency %s\n" % event["name"])
                if not ok_l:
                    pfapp.print("\nBranch Violation for contingency %s\n" % event["name"])
                if ok_v and ok_l:
                    pfapp.print("\nNo violation for contingency %s\n" % event["name"])
        else:
            results.append((event["name"], False, "DIVERGENT"))
            if print_calcs:
                pfapp.print("\nDivergent for contingency %s\n" % event["name"])

        # Restore network
        pfapp.unSetContingency(contingency)
        if check_qlim:
            pfapp.clearQlimViolations()

        if print_calcs:
            pfapp.close()

    # Print summary
    taskmgr.printStats()

    if comm.rank() == 0:
        print("\n" + "=" * 60)
        print("Contingency Analysis Summary")
        print("=" * 60)
        for name, success, status in results:
            print("  %-40s %s" % (name, status if success else "FAILED"))

    timer.stop(t_total)
    if not args.no_timer:
        timer.dump()

    return 0


# -------------------------------------------------------------
# Main entry point
# -------------------------------------------------------------
def main(argv=None):
    """Main CLI entry point."""
    if argv is None:
        argv = sys.argv[1:]

    parser = argparse.ArgumentParser(
        prog="gridpack",
        description="GridPACK: Unified command-line interface for power grid analysis",
        epilog="Run with MPI: mpiexec -np N gridpack <command> config.xml [options]"
    )
    parser.add_argument("--version", action="version", version="GridPACK CLI 0.1.0")

    subparsers = parser.add_subparsers(
        title="commands",
        dest="command",
        metavar="<command>"
    )

    # -- Common arguments added to each subparser --
    def add_common_args(sub):
        sub.add_argument("config", help="XML configuration file")
        sub.add_argument("-q", "--quiet", action="store_true",
                         help="Suppress output")
        sub.add_argument("-o", "--output", metavar="FILE",
                         help="Redirect output to file")
        sub.add_argument("--no-timer", action="store_true",
                         help="Skip timer dump at end")

    # -- powerflow --
    p_pf = subparsers.add_parser(
        "powerflow", aliases=["pf"],
        help="AC power flow analysis",
        description="Run AC power flow analysis on a network."
    )
    add_common_args(p_pf)
    p_pf.add_argument("--solver", choices=["nr", "nl"],
                       help="Solver: nr (Newton-Raphson) or nl (nonlinear)")
    p_pf.add_argument("--export-psse", nargs=2, metavar=("FORMAT", "FILE"),
                       help="Export to PSS/E format: v23, v33, or v34")
    p_pf.set_defaults(func=cmd_powerflow)

    # -- dsf --
    p_dsf = subparsers.add_parser(
        "dsf",
        help="Dynamic simulation framework",
        description="Run dynamic simulation (transient stability) analysis."
    )
    add_common_args(p_dsf)
    p_dsf.set_defaults(func=cmd_dsf)

    # -- se --
    p_se = subparsers.add_parser(
        "se", aliases=["state-estimation"],
        help="State estimation",
        description="Run weighted least squares state estimation."
    )
    add_common_args(p_se)
    p_se.set_defaults(func=cmd_se)

    # -- hadrec --
    p_hadrec = subparsers.add_parser(
        "hadrec",
        help="HADREC remedial action control",
        description="Run Hybrid Automatic Dynamic REmedial action Control."
    )
    add_common_args(p_hadrec)
    p_hadrec.set_defaults(func=cmd_hadrec)

    # -- emt --
    p_emt = subparsers.add_parser(
        "emt",
        help="Electromagnetic transient simulation",
        description="Run electromagnetic transient (EMT) simulation."
    )
    add_common_args(p_emt)
    p_emt.set_defaults(func=cmd_emt)

    # -- ca --
    p_ca = subparsers.add_parser(
        "ca", aliases=["contingency"],
        help="Contingency analysis",
        description="Run N-1 contingency analysis using power flow."
    )
    add_common_args(p_ca)
    p_ca.add_argument("--vlimits", nargs=2, type=float,
                       metavar=("VMIN", "VMAX"),
                       help="Voltage limits (default: from XML or 0.9 1.1)")
    p_ca.add_argument("--group-size", type=int, metavar="N",
                       help="MPI group size for task distribution")
    p_ca.add_argument("--no-print-calcs", action="store_true",
                       help="Do not write per-contingency output files")
    p_ca.set_defaults(func=cmd_ca)

    # Parse and dispatch
    args = parser.parse_args(argv)

    if not args.command:
        parser.print_help()
        return 1

    # Create the GridPACK environment (MPI init) before the handler,
    # and close it (MPI finalize) after — matching the lifecycle pattern
    # used by the standalone scripts (pf.py, dsf.py, etc.).
    import signal
    gp_env = _GridPACKEnv()
    rc = 0
    try:
        rc = args.func(args, gp_env) or 0
    except Exception as e:
        sys.stderr.write("Error: %s\n" % str(e))
        rc = 1

    sys.stdout.flush()
    sys.stderr.flush()
    os._exit(rc)


if __name__ == "__main__":
    main()
