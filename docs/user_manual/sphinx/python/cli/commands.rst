Command Reference
=================

This chapter provides a detailed reference for each CLI subcommand,
including its options, XML configuration, and expected output.

Common Options
--------------

The following options are available across all subcommands:

.. option:: config

   *(Required)* Path to the XML configuration file that defines the
   analysis parameters, network data, and solver settings.

.. option:: -q, --quiet

   Suppress application output. Useful for batch runs where only
   exit status is needed.

.. option:: -o FILE, --output FILE

   Redirect structured output to the specified file. When omitted,
   results are written to standard output.

.. option:: --no-timer

   Omit the timing statistics that are normally printed at the end
   of each run.


.. _cmd-powerflow:

Power Flow Analysis
-------------------

Solve the AC power flow equations for a power system network.

**Subcommand:** ``powerflow`` (alias: ``pf``)

**Synopsis**::

    gridpack powerflow <config.xml> [options]
    gridpack pf <config.xml> [options]

**Subcommand Options:**

.. option:: --solver {nr,nl}

   Override the solver method:

   - ``nr`` — Custom Newton-Raphson iterative solver (default)
   - ``nl`` — PETSc nonlinear solver (SNES)

   When omitted, the solver is selected based on the ``UseNonLinear``
   parameter in the XML configuration.

.. option:: --export-psse FORMAT FILE

   Export the solved network to PSS/E RAW format.
   ``FORMAT`` must be one of ``v23``, ``v33``, or ``v34``.
   ``FILE`` is the output file path.

**XML Configuration:**

The configuration file must contain a ``<Powerflow>`` block under
``<Configuration>`` with the following key parameters:

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Parameter
     - Description
   * - ``networkConfiguration_v33``
     - Path to the PSS/E RAW network data file
   * - ``maxIteration``
     - Maximum number of Newton-Raphson iterations
   * - ``tolerance``
     - Convergence tolerance for power mismatch
   * - ``UseNonLinear``
     - Select nonlinear solver (``true`` / ``false``)
   * - ``exportPSSE_v23`` / ``v33`` / ``v34``
     - Auto-export output filenames
   * - ``suppressOutput``
     - Suppress all output (``true`` / ``false``)

**Output:**

Bus voltage magnitudes and angles, branch power flows, and solver
convergence information.

**Examples**::

    gridpack powerflow input.xml
    gridpack pf input.xml --solver nl --output results.txt
    gridpack pf input.xml --export-psse v33 solved.raw
    mpiexec -np 4 gridpack pf large_network.xml


.. _cmd-dsf:

Dynamic Simulation
------------------

Perform transient stability analysis using the dynamic simulation
framework (DSF).

**Subcommand:** ``dsf``

**Synopsis**::

    gridpack dsf <config.xml> [options]

**XML Configuration:**

The configuration file must contain both ``<Powerflow>`` and
``<Dynamic_simulation>`` blocks under ``<Configuration>``.

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Parameter
     - Description
   * - Network data file
     - PSS/E RAW format (``.raw``)
   * - Dynamic model file
     - PSS/E DYR format (``.dyr``)
   * - Simulation timing
     - Start time, end time, integration time step
   * - Fault definitions
     - Bus faults, line faults, and generator trip events
   * - Generator watch
     - Specifies which generators to monitor during simulation

**Output:**

Time-domain trajectories of generator rotor angles, rotor speeds,
and electrical quantities. Generator watch files are created as
specified in the XML configuration.

**Examples**::

    gridpack dsf input.xml
    gridpack dsf input.xml --quiet --output dsf_results.txt
    mpiexec -np 8 gridpack dsf input.xml


.. _cmd-se:

State Estimation
----------------

Solve the weighted least squares (WLS) state estimation problem.

**Subcommand:** ``se`` (alias: ``state-estimation``)

**Synopsis**::

    gridpack se <config.xml> [options]
    gridpack state-estimation <config.xml> [options]

**XML Configuration:**

The configuration file must contain a ``<State_estimation>`` block
with the following key parameters:

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Parameter
     - Description
   * - ``measurementList``
     - Path to the measurement data XML file
   * - Network data file
     - PSS/E RAW format (``.raw``)

**Measurement Types:**

The measurement XML file supports the following measurement types:

.. list-table::
   :header-rows: 1
   :widths: 15 85

   * - Type
     - Description
   * - ``VM``
     - Bus voltage magnitude
   * - ``VA``
     - Bus voltage angle
   * - ``PIJ``
     - Real power flow on branch (from-bus side)
   * - ``QIJ``
     - Reactive power flow on branch (from-bus side)
   * - ``PJI``
     - Real power flow on branch (to-bus side)
   * - ``QJI``
     - Reactive power flow on branch (to-bus side)
   * - ``IIJ``
     - Current magnitude on branch (from-bus side)
   * - ``IJI``
     - Current magnitude on branch (to-bus side)
   * - ``VUL``
     - Upper voltage limit (pseudo-measurement)
   * - ``VLL``
     - Lower voltage limit (pseudo-measurement)

**Output:**

Estimated bus voltage magnitudes and angles for all buses in the
network. If the estimator does not fully converge, a warning is
printed to standard error; partial results are still written.

**Examples**::

    gridpack se input.xml
    gridpack se input.xml --output se_results.txt


.. _cmd-hadrec:

HADREC
------

Run the Hybrid Automatic Dynamic Remedial Action Control (HADREC)
simulation framework.

**Subcommand:** ``hadrec``

**Synopsis**::

    gridpack hadrec <config.xml> [options]

**Description:**

HADREC integrates power flow initialization and dynamic simulation
with a remedial action control loop. During the simulation, control
actions — including load shedding, line tripping, and generator
tripping — can be applied based on real-time observations of system
state variables such as bus voltages, generator rotor speeds, and
power flows.

**XML Configuration:**

The configuration file should include power flow settings, dynamic
simulation parameters, observation point definitions, and remedial
action specifications.

**Examples**::

    gridpack hadrec input.xml
    gridpack hadrec input.xml --quiet
    mpiexec -np 4 gridpack hadrec input.xml


.. _cmd-emt:

Electromagnetic Transient Simulation
-------------------------------------

Perform electromagnetic transient (EMT) simulation.

**Subcommand:** ``emt``

**Synopsis**::

    gridpack emt <config.xml> [options]

**Description:**

EMT simulation models power system dynamics at the sub-cycle time
scale, resolving fast switching transients, detailed waveforms, and
high-frequency phenomena that are not captured by phasor-domain
(RMS) transient stability analysis.

**Examples**::

    gridpack emt input.xml
    mpiexec -np 4 gridpack emt input.xml


.. _cmd-ca:

Contingency Analysis
--------------------

Perform N-1 contingency analysis by systematically tripping lines
and generators and solving power flow for each contingency case.

**Subcommand:** ``ca`` (alias: ``contingency``)

**Synopsis**::

    gridpack ca <config.xml> [options]
    gridpack contingency <config.xml> [options]

**Subcommand Options:**

.. option:: --vlimits VMIN VMAX

   Override the voltage violation thresholds. When omitted, values
   are read from the XML configuration (default: 0.9 and 1.1 per unit).

.. option:: --no-print-calcs

   Suppress individual per-contingency output files. Only the summary
   table and timing statistics are printed.

**XML Configuration:**

The configuration file must contain two blocks:

1. A ``<Contingency_analysis>`` block with analysis parameters:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Parameter
     - Description
   * - ``contingencyList``
     - Path to the contingency definitions XML file
   * - ``minVoltage``
     - Lower voltage violation threshold (per unit)
   * - ``maxVoltage``
     - Upper voltage violation threshold (per unit)
   * - ``printCalcFiles``
     - Write per-contingency output files (``true`` / ``false``)
   * - ``qlim``
     - Enable reactive power limit checking (``true`` / ``false``)

2. A ``<Powerflow>`` block with network and solver parameters (same
   format as for the power flow subcommand).

**Contingency List Format:**

Contingencies are defined in a separate XML file referenced by the
``contingencyList`` parameter. The file supports two contingency types:

*Line contingency* — trips one or more transmission lines::

    <Contingency>
      <contingencyType>Line</contingencyType>
      <contingencyName>Line_1_2_BL</contingencyName>
      <contingencyLineBuses>1 2</contingencyLineBuses>
      <contingencyLineNames>BL</contingencyLineNames>
    </Contingency>

*Generator contingency* — trips one or more generators::

    <Contingency>
      <contingencyType>Generator</contingencyType>
      <contingencyName>Gen_3_1</contingencyName>
      <contingencyBuses>3</contingencyBuses>
      <contingencyGenerators>1</contingencyGenerators>
    </Contingency>

These elements are wrapped in the following XML structure::

    <ContingencyList>
      <Contingency_analysis>
        <Contingencies>
          <!-- Contingency elements here -->
        </Contingencies>
      </Contingency_analysis>
    </ContingencyList>

**Analysis Workflow:**

The contingency analysis proceeds as follows:

1. Solve the base-case power flow.
2. Record any pre-existing voltage violations so they are excluded
   from contingency violation checks.
3. For each contingency (distributed across MPI processes via the
   TaskManager):

   a. Reset all bus voltages to the base-case solution.
   b. Apply the contingency (trip lines or generators).
   c. Solve power flow for the modified network.
   d. Check for voltage violations and branch overload violations.
   e. Restore the network to its pre-contingency state.

4. Aggregate and print the results summary.

**Output:**

- Per-contingency ``.out`` files containing detailed power flow
  results (unless ``--no-print-calcs`` is specified).
- A summary table reporting convergence status and violation types
  for each contingency case.
- Timing statistics.

**Examples**::

    gridpack ca input.xml
    gridpack ca input.xml --vlimits 0.95 1.05 --no-print-calcs
    mpiexec -np 4 gridpack ca input.xml
