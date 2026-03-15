Overview
========

The GridPACK command-line interface (CLI) provides a unified entry point
for running GridPACK power system analysis applications from the terminal.
It replaces the need to invoke individual Python scripts by consolidating
all supported applications under a single ``gridpack`` command with
intuitive subcommands.

The CLI is implemented in Python and built on the GridPACK pybind11
bindings, providing direct access to the same high-performance C++
solvers used by the native GridPACK executables.

Supported Applications
----------------------

The following table summarizes the available analysis modules. Each
module corresponds to a CLI subcommand and may also have a short alias.

.. list-table::
   :header-rows: 1
   :widths: 20 15 65

   * - Subcommand
     - Alias
     - Description
   * - ``powerflow``
     - ``pf``
     - AC power flow analysis using Newton-Raphson or nonlinear solver
   * - ``dsf``
     -
     - Dynamic simulation framework for transient stability analysis
   * - ``se``
     - ``state-estimation``
     - Weighted least squares state estimation
   * - ``hadrec``
     -
     - Hybrid Automatic Dynamic Remedial Action Control (HADREC)
   * - ``emt``
     -
     - Electromagnetic transient simulation
   * - ``ca``
     - ``contingency``
     - N-1 contingency analysis with parallel task distribution

Command Structure
-----------------

All commands follow a consistent invocation pattern::

    gridpack <subcommand> <config.xml> [options]

For parallel execution using MPI::

    mpiexec -np <N> gridpack <subcommand> <config.xml> [options]

Each subcommand accepts the XML configuration file as its first positional
argument. The XML file specifies all application-specific parameters, including
network data paths, solver settings, and output options. Additional command-line
flags allow overriding selected parameters without editing the XML.

Relationship to Existing Scripts
--------------------------------

The CLI does not replace the standalone Python scripts (``pf.py``,
``dsf.py``, ``stes.py``, ``hadrec.py``, ``emt.py``) that ship with
GridPACK. These scripts remain available and continue to work as before.
The CLI provides an additional, unified interface that is better suited
for production use, shell scripting, and Docker-based workflows.

Getting Help
------------

To display the list of available subcommands::

    gridpack --help

To display options for a specific subcommand::

    gridpack powerflow --help
    gridpack ca --help
