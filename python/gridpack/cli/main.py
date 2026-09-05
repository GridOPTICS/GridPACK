#!/usr/bin/env python
# -------------------------------------------------------------
# file: gridpack/cli/main.py
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

Exit codes:
    0  success
    1  unexpected error
    2  bad arguments or unreadable config
    3  the analysis did not converge

A thin argparse shell over the ``gridpack.*`` wrappers.  Reaching into
``gridpack._gridpack`` here duplicates library logic (that is how the
XML-bool and divergence bugs got in), so ``tests/test_cli_uses_library.py``
fails the build on new low-level calls; CoarseTimer and emt.EMT are the
documented exceptions.
"""

from __future__ import print_function
import sys
import os
import argparse


EXIT_OK = 0
EXIT_ERROR = 1
EXIT_USAGE = 2
EXIT_DIVERGED = 3

_PSSE_FORMATS = {"v23": 23, "v33": 33, "v34": 34}


# -------------------------------------------------------------
# Powerflow subcommand
# -------------------------------------------------------------
def cmd_powerflow(args, session):
    """Run AC power flow analysis."""
    import gridpack

    if args.export_psse:
        fmt = args.export_psse[0]
        if fmt not in _PSSE_FORMATS:
            sys.stderr.write(
                "Error: --export-psse format must be one of %s, not %r\n"
                % (", ".join(sorted(_PSSE_FORMATS)), fmt))
            return EXIT_USAGE

    timer = gridpack.CoarseTimer()

    # None leaves the XML in charge; --solver overrides UseNonLinear.
    nonlinear = {"nl": True, "nr": False}.get(args.solver)

    pf = gridpack.PowerFlow(session, args.config,
                            suppress_output=True if args.quiet else None)

    # strict=False: the divergence verdict is this command's exit code, and
    # the results are still written either way.
    result = pf.solve(nonlinear=nonlinear, strict=False)
    result.write(args.output)

    if args.export_psse:
        fmt, fname = args.export_psse
        result.export_psse(fname, version=_PSSE_FORMATS[fmt])
    for version, path in pf.psse_exports_from_xml:
        result.export_psse(path, version=version)

    if not pf.suppress_output and not args.no_timer:
        timer.dump()

    if not result.converged:
        sys.stderr.write(
            "Error: power flow did not converge; results were written anyway.\n")
        return EXIT_DIVERGED
    return EXIT_OK


# -------------------------------------------------------------
# Dynamic Simulation subcommand
# -------------------------------------------------------------
def cmd_dsf(args, session):
    """Run dynamic simulation framework."""
    import gridpack
    from gridpack.dynamic_simulation import DSFullApp

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
def cmd_se(args, session):
    """Run state estimation."""
    import gridpack

    # SEApp has no open()/close(), so -o used to raise AttributeError
    # after the solve had already printed to stdout.
    if args.output:
        sys.stderr.write(
            "Error: se cannot redirect output to a file; the SEApp binding "
            "exposes no open()/close(). Redirect stdout instead.\n")
        return EXIT_USAGE

    timer = gridpack.CoarseTimer()

    se = gridpack.StateEstimation(session, args.config,
                                  suppress_output=args.quiet)
    result = se.solve()
    result.write()

    if not args.no_timer:
        timer.dump()

    if not result.has_converged():
        sys.stderr.write("Error: state estimation did not converge; "
                         "results were written anyway.\n")
        return EXIT_DIVERGED
    return EXIT_OK


# -------------------------------------------------------------
# HADREC subcommand
# -------------------------------------------------------------
def cmd_hadrec(args, session):
    """Run HADREC remedial action control simulation."""
    import gridpack
    import gridpack.hadrec
    import gridpack.dynamic_simulation

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
def cmd_emt(args, session):
    """Run electromagnetic transient simulation."""
    import gridpack
    import gridpack.emt

    emt_app = gridpack.emt.EMT()
    emt_app.setconfigurationfile(args.config)
    emt_app.solvepowerflow()
    emt_app.setup()
    emt_app.solve()

    return 0


# -------------------------------------------------------------
# Contingency Analysis subcommand
# -------------------------------------------------------------
def cmd_ca(args, session):
    """Run N-1/N-k contingency analysis."""
    import gridpack

    timer = gridpack.CoarseTimer()
    t_total = timer.createCategory("CA: Total Application")
    timer.start(t_total)

    ca = gridpack.ContingencyAnalysis(
        session, args.config,
        voltage_limits=args.vlimits,
        print_calc_files=False if args.no_print_calcs else None,
        suppress_output=True if args.quiet else None,
    )

    # Analyzing nothing must not look like a clean run.  Before the base
    # case, so a config error costs no solve.
    if not ca.contingencies:
        if ca.requests_auto_n1:
            sys.stderr.write(
                "Error: FullBranchN1/FullGeneratorN1 auto-generation is not "
                "supported; list contingencies explicitly via "
                "<contingencyList>.\n")
        else:
            sys.stderr.write(
                "Error: no contingencies to run; set <contingencyList> in the "
                "Contingency_analysis block.\n")
        return EXIT_USAGE

    if session.rank == 0:
        print("=" * 60)
        print("Contingency Analysis")
        print("Total contingencies: %d" % len(ca.contingencies))
        print("Voltage limits: [%.3f, %.3f]" % (ca.min_voltage, ca.max_voltage))
        print("=" * 60)

    # Contingencies are measured against the base case, so a diverged base
    # case makes all of them meaningless.
    if not ca.solve_base_case().converged:
        sys.stderr.write("Error: contingency analysis base case did not "
                         "converge; no contingencies were run.\n")
        return EXIT_DIVERGED

    def announce(task_id, contingency):
        if session.rank == 0:
            print("Running contingency %d: %s" % (task_id, contingency.name))

    ca.run(progress=announce)

    # Every rank evaluated a disjoint subset; gather for a complete summary.
    if session.rank == 0:
        print("\n" + "=" * 60)
        print("Contingency Analysis Summary")
        print("=" * 60)
    results = ca.gather()
    if session.rank == 0:
        for r in results:
            print("  %-40s %s" % (r.name, r.status))

    timer.stop(t_total)
    if not args.no_timer:
        timer.dump()

    # A contingency that diverges is a finding, not a run failure.
    return EXIT_OK


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
    p_ca.add_argument("--no-print-calcs", action="store_true",
                       help="Do not write per-contingency output files")
    p_ca.set_defaults(func=cmd_ca)

    # Parse and dispatch
    args = parser.parse_args(argv)

    if not args.command:
        parser.print_help()
        return EXIT_USAGE

    # Checked here so a typo exits cleanly instead of failing inside the
    # C++ config reader after MPI is already up.
    if not os.path.isfile(args.config):
        sys.stderr.write("Error: no such config file: %s\n" % args.config)
        return EXIT_USAGE

    import gridpack

    session = gridpack.Session()
    rc = EXIT_OK
    try:
        rc = args.func(args, session) or EXIT_OK
    except Exception as e:
        # GRIDPACK_TRACEBACK=1 keeps the trace, which is otherwise the
        # only way to diagnose a failure on one rank.
        if os.environ.get("GRIDPACK_TRACEBACK"):
            import traceback
            traceback.print_exc()
        sys.stderr.write("Error: %s: %s\n" % (type(e).__name__, e))
        rc = EXIT_ERROR
    finally:
        sys.stdout.flush()
        sys.stderr.flush()

    # Session.close() drains its registered wrappers before releasing the
    # Communicator and Environment, which is what lets MPI finalize.  This
    # used to be os._exit(rc) to dodge a DSFullApp teardown SEGV, but that
    # skipped MPI_Finalize, so mpiexec reported abnormal termination (exit 1)
    # even on a clean run.  The SEGV no longer reproduces; the exit code did.
    session.close()
    return rc


if __name__ == "__main__":
    sys.exit(main())
