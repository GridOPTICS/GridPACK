Application Modules
===================

Many of the example applications in GridPACK have been converted to
modules that can be called from other programs. These modules make it
relatively simple to chain different types of calculations together to
form larger applications. An example is using power flow or state
estimation to initialize a dynamic simulation. The modules are designed
to separat the major phases of the calculation into calls so that users
have some fine-grained control that allows them to mix different
applications together. In most cases, different options for setting up
calculations are provided so that once a network has been read in and
partitioned, it is not necessary to repeat this process when a new
calculation is started based on the results of a previous simulation.

Currently, four applications are available as modules within GridPACK.
They include power flow, state estimation, dynamic simulation using
the full Y-matrix, and electromagnetic transient (EMT) simulation. Each of 
these modules can be used to create a short, standalone application, but the 
goal is to enable users to combine modules together in more complicated work 
flows. These modules can also be used as a starting point for users to create 
their own applications by modifying the existing code in the modules to create 
new functionality. Each of the modules is described in more detail below.
Example codes that use the modules to implement applications can be
found in the ``src/application`` directory. These include powerflow,
state estimation, contingency analysis, dynamic simulation, and EMT simulation. 
These directories also contain sample input networks and input files. Options
for different solvers can be found in these files.

Most of these applications can use different format PSS/E files for input. The
standard format is PSS/E v23 and is specified via the ``networkConfiguration``
field in the input XML file. The v33, v34, v35 and v36 formats can be used by
appending ``_v33``, ``_v34``, ``_v35``, or ``v_36`` to the
``networkConfiguration`` field. For example, to include a version 33
network configuration
file, modify the network configuration line in the input XML file to

::

   <networkConfiguration_v33> v33_format.raw </networkConfiguration_v33>

Power Flow
----------

The power flow module consists of a collection of function calls that
can be used to set up and run power flow calculations. Additional
routines are designed to support different types of contingency
analysis. The power flow application class is ``PFAppModule`` and
belongs to the ``gridpack::powerflow`` namespace. The constructor
and destructor for this class are simple and only create the basis power
flow object. In particular, the power flow network must be created
outside the power flow object and then assigned to the object when the
network configuration file is read in. This can be done with the call

::

   void readNetwork(boost::shared_ptr<PFNetwork> &network,
         Configuration *config)

The ``Configuation`` object should already be pointing to an open
file containing a ``Powerflow`` block. This block contains a
``networkConfiguration`` field that has the name of the PSS/E format
file containing the network information. An example powerflow input file is
shown below:

::

  <?xml version="1.0" encoding="utf-8"?>
  <Configuration>
    <Powerflow>
      <networkConfiguration> IEEE14.raw </networkConfiguration>
      <maxIteration>50</maxIteration>
      <tolerance>1.0e-6</tolerance>
      <dampingFactor>1.0</dampingFactor>
      <initStart>warm</initStart>
      <qlim>true</qlim>
      <maxQlimIterations>3</maxQlimIterations>
      <qlimDeadband>0.1</qlimDeadband>
      <SwitchedShunt>false</SwitchedShunt>
      <LTC>false</LTC>
      <AreaInterchange>false</AreaInterchange>
      <maxControllerIterations>10</maxControllerIterations>
      <outputFormat>json</outputFormat>
      <outputFile>pf_result</outputFile>
      <LinearSolver>
        <PETScOptions>
          -ksp_type preonly
          -pc_type lu
          -pc_factor_mat_solver_type superlu_dist
        </PETScOptions>
      </LinearSolver>
    </Powerflow>
  </Configuration>

The ``networkConfiguration`` field accepts PSS/E RAW files in any
supported version (v23, v33, v34, v35, v36). For v30+ files, the parser
automatically detects the version from the third field of the first line
in the RAW file header. Version-specific tags
(``networkConfiguration_v33``, ``networkConfiguration_v34``, etc.) are
still supported for backward compatibility and take precedence over
auto-detection. PSS/E v23 files (which lack the version field) are
handled as the default when auto-detection finds no version indicator.

This example specifies the input network configuration, the maximum number of
iterations in the non-linear Newton-Raphson solver, the solution tolerance
and the properties of the linear solver. The ``qlim`` parameter enables
reactive power limit enforcement (default ``true``). When enabled, PV
buses whose generator reactive power output exceeds the specified
``Qmax`` or ``Qmin`` limits are switched to PQ buses and the power
flow is re-solved. The ``maxQlimIterations`` parameter (default 3)
controls the maximum number of outer PV-to-PQ switching iterations.
The ``qlimDeadband`` parameter (default 0.1 Mvar) sets a tolerance around
the reactive power limits: a PV bus is only converted to PQ when its
required Q exceeds a limit by more than this amount, preventing spurious
switching due to small numerical imbalances near the limit.

The ``SwitchedShunt`` parameter (default ``false``) enables automatic
switched shunt voltage control. When enabled, buses with switched shunt
data (MODSW=1 for discrete or MODSW=2 for continuous control) are
monitored after each Newton-Raphson convergence. If the controlled bus
voltage falls outside the ``[VSWLO, VSWHI]`` deadband, the shunt
susceptance is adjusted: discrete mode switches one capacitor or reactor
bank step per iteration, while continuous mode adjusts B smoothly toward
the deadband midpoint. Cycle detection locks shunts that oscillate
between two states. Shunts with ``SWREM`` set to a nonzero bus number
regulate voltage at that remote bus instead of the local bus.

The ``LTC`` parameter (default ``false``) enables load tap changer
(LTC) control on transformers. Transformers with COD1=1 in their PSS/E
data (or a nonzero CONT bus in v23 transformer adjustment records)
automatically adjust their tap ratios to regulate voltage at the
controlled bus within the ``[VMI, VMA]`` deadband. Tap ratios are
adjusted by one discrete step per controller iteration, bounded by
``[RMI, RMA]``. The step size is computed from the tap range and
number of tap positions (NTP), or read directly from the STEP field
when available. Cycle detection prevents tap hunting.

Remote voltage regulation (IREG) is enabled automatically for any
generator whose PSS/E ``IREG`` field points to a bus other than its
own. The controlled remote bus is treated as PV (held at the
generator's scheduled voltage ``VS``) and the local bus is solved
normally; the source bus voltage is adjusted across the controller loop
until the remote bus tracks setpoint. No XML option is required. When
the remote bus is the swing bus the regulation is dropped with a
warning, and when multiple generators on the same source bus include a
locally-regulating unit (``IREG=0``) local control wins.

The ``AreaInterchange`` parameter (default ``false``) enables area
interchange control. When enabled, after the inner controller loop
converges, the solver computes actual MW exports for each area by
summing tie-line flows (branches crossing area boundaries). If the
export deviates from the desired value (``PDES`` from the area
interchange data) by more than the tolerance (``PTOL``), the solver
adjusts real power generation at the area slack bus (``ISW``) and
re-solves. This outer loop runs up to 10 iterations.

The ``maxControllerIterations`` parameter (default 10) sets the maximum
number of inner controller iterations. This limit is shared by all
active controls (Q-limits, switched shunts, LTC, and IREG remote
voltage regulation). When only Q-limits are active, the
``maxQlimIterations`` parameter is used for backward compatibility.

The ``dampingFactor`` parameter (default 1.0) scales each Newton-Raphson
correction vector before it is applied to the bus voltages and angles.
A value of 1.0 gives standard (undamped) Newton-Raphson. Values between
0 and 1 reduce the step size, which can improve convergence stability for
ill-conditioned systems such as heavily loaded networks or cases where the
warm-start initialization is far from the solution. For well-initialized
cases (RAW files with converged VM/VA stored) a value of 1.0 is
recommended.

The ``initStart`` parameter (default ``warm``) controls how bus voltages
and angles are initialized before the first Newton-Raphson iteration.
``warm`` reads voltage magnitude and angle from the PSS/E RAW file bus
section (columns VM and VA), placing all buses at the operating point
stored in the network file. ``flat`` initializes all buses to 1.0 pu and
0 degrees regardless of the RAW file values, which can be useful for
cases where the stored voltages are unreliable or absent.

The ``outputFormat`` and ``outputFile`` parameters control structured
output after a successful solve. When ``outputFormat`` is set to
``json``, the solver writes per-bus voltages, angles, and power
injections to ``<outputFile>_buses.json`` and per-branch flows to
``<outputFile>_branches.json``. When set to ``csv``, the results are
written to ``<outputFile>_buses.csv`` and ``<outputFile>_branches.csv``.
If these parameters are omitted no structured output file is written.

The recommended ``LinearSolver`` configuration for the hand-coded
Newton-Raphson ``solve()`` routine is ``-ksp_type preonly`` with direct
LU factorization (``-pc_type lu``). This applies the factorization once
per iteration without an iterative Krylov loop and avoids spurious
divergence checks that can terminate otherwise-converging solutions.

This is a minimal example and more options for the solver can be used.
Readers are encouraged to look at the examples that are copied into the
powerflow directory as part of the build.

The network configuration file
is read directly from the input deck by the
``readNetwork`` method. The ``PFNetwork`` is defined in
the the ``gridpack.hpp`` header file. The configuration module is
usually opened in the main calling program and a pointer to the file can
be passed through to power flow module. The ``readNetwork`` routine
also partitions the network.

Once the network has been read in, the internal indices and exchange
buffers can be set up by calling

::

   void initialize()

The power flow application is now ready to be used. To solve the current
configuration, the calls

::

   void solve()

   void nl_solve()

can be used. The first call solves the system uses a hand coded
Newton-Raphson iteration loop to solve the system, the second call uses
a non-linear solver to solve the power flow equations. Both solvers can
be controlled through solver options in the input file. The type of
linear solver used in the solve routine is controlled by the parameters
in the ``LinearSolver`` block, the non-linear solver is controlled
by the properties in the ``NonlinearSolver`` block Output from the
power flow solution can be written to an output file or standard out
using one of the commands

::

   void write()

   void writeBus(const char* signal)

   void writeBranch(const char* signal)

The first command writes out the real and imaginary parts of the complex
power for the branches and the voltage magnitude and phase angle for the
buses. The second command only writes out bus properties. If no argument
is given, the command writes out the voltage magnitude and phase angle
for every bus. For buses, the argument ``pq`` writes out the real
and imaginary parts of the complex voltage and ``record`` writes
out the type of bus, the total active and reactive constant power loads,
and the total active and reactive generator power outputs. For branches,
``flow`` writes out the real and imaginary parts of the complex
power and “``record``” writes out the values of the resistance,
reactance, charging and A, B, C ratings for each line element.

Additional information can be written to standard out or a file using
the command

::

   void print(const char* buf)

which writes out the contents of the character array ``buf``. This
command can be called from all processors, but only one processor
actually writes out data.

The location of output can be controlled using the commands

::

   void open(const char* filename)

   void close()

If the write commands or ``print`` are used without calling
``open``, then all output is directed to standard out. If
``open`` is called, then the output is directed to the file
specified in filename until the ``close`` command is called, after
which all output is again directed towards standard out.

If the results of the power flow calculation are needed by another
calculation, then the voltage magnitude and phase angle of the bus and
the real and imaginary parts of the complex power for each generator can
be stored in the ``DataCollection`` objects on each bus using the
command

::

   void saveData()

If the network is then copied to a new type of network, this information
is carried over to the new network. The voltage magnitude and phase
angle is stored in the ``DataCollection`` variables
``BUS_PF_VMAG`` and ``BUS_PF_VANG`` and the generator parameters
are stored in the indexed variables ``GENERATOR_PF_PGEN[i]`` and
``GENERATOR_PF_QGEN[i]``, where the index ``i`` runs over all
generators on the bus.

The power flow module supports ZIP loads, where the total load power is
modeled as a combination of constant power (P), constant current (I),
and constant admittance (Y) components. The load power at voltage
``V`` is computed as ``P_load = PL + IP*V + YP*V^2`` and
``Q_load = QL + IQ*V - YQ*V^2``, following PSS/E sign conventions.
The ZIP load coefficients are read automatically from PSS/E network
files when present.

Switched shunt devices specified in PSS/E input files are supported.
The initial susceptance value (``BINIT``) is absorbed into the fixed
bus shunt admittance during network setup, consistent with the behavior
of commercial power system tools.

The power flow module includes island detection capabilities. The
function ``getIslandCount()`` returns the number of electrical
islands found in the network after calling ``solve()``. The function
``hasLoneBus()`` identifies isolated buses with no active connections.
The function ``checkAndTransferSlack()`` automatically transfers the
slack bus role to the bus with the largest online generator capacity if
the original slack bus generator goes offline.

The output from ``write()`` includes per-bus power injection values
(``Pinj`` and ``Qinj``) in addition to voltage magnitude and phase
angle. These values reflect the net power injection at each bus
including all ZIP load contributions at the solved voltage.

The remaining methods in the PFAppModule class support different kinds
of contingency applications. Contingencies are defined using the data
structure

::

   struct Contingency
   {
     int p_type;
     std::string p_name;
     // Line contingencies
     std::vector<int> p_from;
     std::vector<int> p_to;
     std::vector<std::string> p_ckt;
     // Status of line before contingency
     std::vector<bool> p_saveLineStatus;
     // Generator contingencies
     std::vector<int> p_busid;
     std::vector<std::string> p_genid;
     // Status of generator before contingency
     std::vector<bool> p_saveGenStatus;
   };

The variable ``p_type`` corresponds to an enumerated type that can
have the values ``Generator`` and ``Branch``. The variables
``p_saveLinesStatus`` and ``p_saveGenStatus`` are used
internally and do not have to be set by the user. The remaining
variables are used to describe the lines and generators that may fail
during a contingency event. These variables are all vectors, since a
single contingency could theoretically represent the failure of multiple
elements. For failures of type ``Branch``, the variables
``p_from`` and ``p_to`` are the original indices of the “from”
and “to” bus that identify a branch and the variable ``p_ckt`` is
the 2 character identifier of the individual transmission element. For
failures of type ``Generator``, ``p_busid`` is the original
index of the bus and ``p_genid`` is the 2 character identifier of
the generator that fails. An example of how to use this functionality is
given in the contingency analysis application that can be found under
``src/applications/contingency_analysis``. This is also a good
example of how to use modules.

Two calls

::

   bool setContingency(Contingency &event)

   bool unSetContingency(Contingency &event)

can be used to set or unset a contingency. Both functions return
``true`` if the contingency was successfully applied or removed. The
call ``unSetContingency`` should only be called after calling
``setContingency`` and it should use the same ``event``
argument. After calling the ``unSetContingency`` method, the network
should have the same configuration as before calling the
``setContingency`` method.

After applying a contingency, the network topology should be checked
before attempting to solve. The function

::

   int getIslandCount()

returns the number of electrical islands in the network. A value greater
than 1 indicates that the contingency has split the network and the
power flow solve should be skipped. The function

::

   bool hasLoneBus()

returns ``true`` if any bus has become completely isolated with no
active branch connections.

If a contingency takes the slack bus generator offline, the slack
designation can be transferred using

::

   bool checkAndTransferSlack()

This function transfers the slack bus role to the bus with the largest
online generator capacity. It returns ``false`` if no valid slack bus
can be found, indicating the system is unsolvable. After unSetting a
contingency, the original slack bus can be restored by calling

::

   void restoreSlack()

After solving the post-contingency power flow, the slack bus generation
can be validated using

::

   bool checkSlackCapacity()

This returns ``false`` if the slack bus generator output exceeds its
``Pmax`` limit, indicating insufficient generation capacity.

The remaining calls in ``PFAppModule`` can be used to determine the
status of a network after solving a configuration with a contingency.
Voltage limits for violation checking are set using

::

   void setVoltageLimits(double Vmin, double Vmax)

where ``Vmin`` and ``Vmax`` are the minimum and maximum allowable
voltage excursions (per unit). The functions

::

   bool checkVoltageViolations()

   bool checkVoltageViolations(int area)

can then be used to check for voltage violations anywhere in the system
or only on buses with the specified value of ``area``. These
functions return ``true`` if there are no voltage violations and
``false`` if a violation is found on one or more buses. It frequently
turns out that many networks have voltage violations even in the absence
of any contingencies and it is often desirable to ignore these
violations. This can be accomplished using the function

::

   void ignoreVoltageViolations()

If this function is called after solving the power flow system in the
absence of any contingencies, then buses that contain violations will be
ignored in subsequent checks of violations. These settings can be undone
by calling

::

   void clearVoltageViolations()

Line overload violations can be checked by calling one of the functions

::

   bool checkLineOverloadViolations()

   bool checkLineOverloadViolations(int area)

   bool checkLineOverloadViolations(std::vector<int> &bus1,
       std::vector<int> &bus2, std::vector<std::string> &tags,
       std::vector<bool> &violations)

The limits on the line are contained in parameters read in from the
network configuration file so the first two functions have no arguments
describing the line limits. The second function only checks for
violations on lines with the specified value of ``area``. The third
overload checks a specific set of branches identified by their endpoint
buses and circuit tags, and populates a per-branch ``violations``
vector. Like voltage violations, branches that display line overload
violations that are present even without contingencies can be ignored in
the checks by calling the function

::

   void ignoreLineOverloadViolations()

after running a calculation on the system without contingencies. These
settings can be cleared using the function

::

   void clearLineOverloadViolations()

Reactive power limit violations on PV buses can be checked using

::

   bool checkQlimViolations()

   bool checkQlimViolations(int area)

When a PV bus generator exceeds its reactive power limits, the bus is
converted to a PQ bus. After unsetting a contingency, these conversions
can be reversed by calling

::

   void clearQlimViolations()

Switched shunt voltage violations can be checked and adjusted using

::

   bool checkSwitchedShuntViolations()

This adjusts shunt susceptance on buses where the controlled voltage is
outside the deadband. After unsetting a contingency, shunt adjustments can
be reset to their initial BINIT values by calling

::

   void clearSwitchedShunts()

LTC transformer tap ratio violations can be checked and adjusted using

::

   bool checkLTCViolations()

This adjusts tap ratios on LTC-controlled transformers where the controlled
bus voltage is outside the deadband. After unsetting a contingency, tap
ratios can be reset to their initial values by calling

::

   void clearLTCControls()

Finally, the internal voltage variables that are used as the solution
variables in the power flow calculation can be reset to their original
values (specified in the network configuration file) by calling the
function

::

   void resetVoltages()

Again, this may be useful in contingency calculations where multiple
calculations are run on the same network and it is desirable that they
all start with the same initial condition.

Contingency Analysis Module
---------------------------

The production contingency analysis driver lives in
``src/applications/contingency_analysis``. It is built on
``PFAppModule`` and adds an MPI task manager, a contingency parser,
monitor filtering, and per-(contingency, branch) CSV output for
downstream statistical analysis. Configuration goes under a
``Contingency_analysis`` block in the input deck alongside the standard
``Powerflow`` block. The XML reader for the contingency list itself is
covered in Section 10; this section documents the runtime options.

Input options
~~~~~~~~~~~~~

+----------------------------+--------+--------------------------------------------------------------+
| Option                     | Default| Description                                                  |
+============================+========+==============================================================+
| ``contingencyList``        | (none) | Path to the contingency XML. May be combined with            |
|                            |        | ``FullBranchN1`` / ``FullGeneratorN1`` for N-1 + custom N-K. |
+----------------------------+--------+--------------------------------------------------------------+
| ``FullBranchN1``           | false  | Auto-generate N-1 over every in-service branch.              |
+----------------------------+--------+--------------------------------------------------------------+
| ``FullGeneratorN1``        | false  | Auto-generate N-1 over every in-service generator.           |
+----------------------------+--------+--------------------------------------------------------------+
| ``groupSize``              | 1      | Deprecated; ignored (forced to 1). An outaged branch can     |
|                            |        | straddle a multi-rank partition. Add ranks for parallelism. |
+----------------------------+--------+--------------------------------------------------------------+
| ``minVoltage``             | 0.9    | Lower voltage limit (pu) for violation checks.               |
+----------------------------+--------+--------------------------------------------------------------+
| ``maxVoltage``             | 1.1    | Upper voltage limit (pu) for violation checks.               |
+----------------------------+--------+--------------------------------------------------------------+
| ``qlim``                   | false  | Enable PV→PQ Q-limit enforcement during the contingency      |
|                            |        | solve. Honors ``qlimDeadband`` from the ``Powerflow`` block. |
+----------------------------+--------+--------------------------------------------------------------+
| ``printCalcFiles``         | true   | Write per-contingency text output (``<name>.out``).          |
+----------------------------+--------+--------------------------------------------------------------+
| ``outputFormat``           | text   | ``text`` | ``json`` | ``csv`` | ``csv_flat`` | ``csv_delta``.|
+----------------------------+--------+--------------------------------------------------------------+
| ``outputFile``             |        | Base name (prefix) for the structured output files.          |
|                            | results|                                                              |
+----------------------------+--------+--------------------------------------------------------------+
| ``writeStats``             | true   | Emit the StatBlock summary ``.txt`` files. Set ``false`` to  |
|                            |        | skip the per-case StatBlock work when CSV is sufficient.     |
+----------------------------+--------+--------------------------------------------------------------+
| ``contingencyRating``      | A      | Loading% denominator across every CA output. ``A``, ``B``,   |
|                            |        | or ``C`` with A→B→C fallback. Default ``A`` matches PW /     |
|                            |        | PSS/E ACCC. ``base_rate_mva`` always uses rate-A.            |
+----------------------------+--------+--------------------------------------------------------------+
| ``monitorBranchesFile``    | (none) | Path to a CSV allowlist (``from_bus,to_bus,ckt`` rows).      |
|                            |        | When set, area/kV gates are ignored.                         |
+----------------------------+--------+--------------------------------------------------------------+
| ``monitorAreas``           | (none) | Area numbers (space- or comma-separated). Branch passes if   |
|                            |        | **either endpoint** is in the set (catches tie-lines).       |
+----------------------------+--------+--------------------------------------------------------------+
| ``monitorKvMin``           | 0      | Lower kV threshold; branch passes if                         |
|                            |        | ``max(kv_from, kv_to) >= monitorKvMin``. 0 disables.         |
+----------------------------+--------+--------------------------------------------------------------+
| ``monitorKvMax``           | 0      | Upper kV threshold; branch passes if                         |
|                            |        | ``max(kv_from, kv_to) <= monitorKvMax``. 0 disables.         |
+----------------------------+--------+--------------------------------------------------------------+

Violation / ranking options driving ``_violations.csv`` and
``_summary.json``:

* ``violationSeverityThreshold`` (``1.0``) — loading% > threshold × 100
  flags a branch violation.
* ``topN`` (``20``, clamped ``[1,10000]``) — cap on ranked arrays.
* ``piBranchWeight`` / ``piVoltageWeight`` (both ``1.0``) — weights on
  the branch-PI and voltage-PI terms in the composite PI ranking.

Monitor filter precedence: When ``monitorBranchesFile`` is set, area/kV 
options are ignored and a warning is logged. Otherwise ``monitorAreas`` 
AND the kV bounds combine.  With no filter set every branch is monitored. 
Circuit IDs are matched with whitespace trimmed on both sides so PSS/E ckt
strings with leading or trailing spaces compare equal to the unpadded form.

A complete annotated example is shipped at
``src/applications/data_sets/input/ca/input_14_filters_example.xml``
with a sample monitor list ``monitor_branches_14.csv``.

CSV outputs (``outputFormat=csv_flat`` / ``csv_delta``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When ``outputFormat`` is ``csv_flat`` or ``csv_delta`` the driver writes
per-(contingency, monitored branch) rows for downstream analysis instead
of the aggregated ``.txt`` files. All file names use the value of
``outputFile`` as a prefix; the default prefix is ``ca_results``. Both
formats also work for generator contingencies — column 3 (``type``)
labels what was tripped and column 31 (``cont_event_facility``)
identifies the tripped element.

``<outputFile>_delta.csv`` (``csv_delta`` only) — wide form, one row
per (contingency, monitored branch) joining base and contingency state
on the same row. Most downstream consumers prefer this because each row
is self-contained.

+-----+----------------------------+----------------------------------------------------------+
| #   | Column                     | Notes                                                    |
+=====+============================+==========================================================+
| 1   | ``event_idx``              | 0 = base case, 1..N = contingencies in input-deck order. |
+-----+----------------------------+----------------------------------------------------------+
| 2   | ``contingency``            | Contingency name (``base_case`` for the base row).       |
+-----+----------------------------+----------------------------------------------------------+
| 3   | ``type``                   | ``branch`` or ``generator`` — what was tripped.          |
+-----+----------------------------+----------------------------------------------------------+
| 4–6 | ``from_bus``, ``to_bus``,  | Identity of the **monitored branch** (not the tripped    |
|     | ``ckt``                    | element).                                                |
+-----+----------------------------+----------------------------------------------------------+
| 7–8 | ``base_kv_from``,          | Endpoint base kV.                                        |
|     | ``base_kv_to``             |                                                          |
+-----+----------------------------+----------------------------------------------------------+
| 9–10| ``area_from``,             | PSS/E area numbers.                                      |
|     | ``area_to``                |                                                          |
+-----+----------------------------+----------------------------------------------------------+
| 11  | ``base_rate_mva``          | Always rate-A (PSS/E "normal" rating).                   |
+-----+----------------------------+----------------------------------------------------------+
| 12  | ``cont_rate_mva``          | Rating selected by ``contingencyRating``.                |
+-----+----------------------------+----------------------------------------------------------+
| 13–14| ``base_p_mw``,            | Real-power flow before / after.                          |
|     | ``cont_p_mw``              |                                                          |
+-----+----------------------------+----------------------------------------------------------+
| 15–16| ``base_q_mvar``,          | Reactive-power flow before / after.                      |
|     | ``cont_q_mvar``            |                                                          |
+-----+----------------------------+----------------------------------------------------------+
| 17–18| ``base_mva``,             | ``sqrt(P² + Q²)`` before / after.                        |
|     | ``cont_mva``               |                                                          |
+-----+----------------------------+----------------------------------------------------------+
| 19  | ``base_loading_pct``       | ``base_mva / base_rate_mva × 100``.                      |
+-----+----------------------------+----------------------------------------------------------+
| 20  | ``cont_loading_pct``       | ``cont_mva / cont_rate_mva × 100``.                      |
+-----+----------------------------+----------------------------------------------------------+
| 21–22| ``v_from_base``,          | From-bus voltage magnitude (pu).                         |
|     | ``v_from_cont``            |                                                          |
+-----+----------------------------+----------------------------------------------------------+
| 23–24| ``v_to_base``,            | To-bus voltage magnitude (pu).                           |
|     | ``v_to_cont``              |                                                          |
+-----+----------------------------+----------------------------------------------------------+
| 25–28| ``ang_from_base``,        | Endpoint bus angles (deg).                               |
|     | ``ang_from_cont``,         |                                                          |
|     | ``ang_to_base``,           |                                                          |
|     | ``ang_to_cont``            |                                                          |
+-----+----------------------------+----------------------------------------------------------+
| 29–30| ``d_angle_base``,         | ``ang_from − ang_to`` before / after.                    |
|     | ``d_angle_cont``           |                                                          |
+-----+----------------------------+----------------------------------------------------------+
| 31  | ``cont_event_facility``    | The tripped element (e.g. ``[area] from to ckt`` for a   |
|     |                            | branch trip, ``gen <bus> <id>`` for a gen trip).         |
+-----+----------------------------+----------------------------------------------------------+

``<outputFile>_flat.csv`` (``csv_flat`` only) — long form, one row per
(case, branch). Columns: ``event_idx, contingency, from_bus, to_bus,
ckt, p_from_mw, q_from_mvar, mva_from, rate_mva, loading_percent, viol,
v_from_pu, v_to_pu, ang_from_deg, ang_to_deg``. ``rate_mva`` is rate-A
on the base row and the configured ``contingencyRating`` on contingency
rows.

``<outputFile>_buses.csv`` (both CSV formats) — bus metadata sidecar so
the per-branch files can stay narrow. Columns: ``bus_id, bus_name,
base_kv, area, zone, owner, area_name, zone_name, owner_name``.

``<outputFile>_convergence.csv`` (every ``outputFormat``) — one row per
contingency. Columns: ``event_idx, contingency, type, converged,
iterations, final_tolerance, max_p_bus, max_p_mismatch, max_q_bus,
max_q_mismatch, status_code``. ``converged`` is ``true`` iff
``status_code == "OK"``; ``status_code`` is one of ``OK`` / ``ISLANDED``
/ ``NO_SLACK`` / ``DIVERGED`` / ``SLACK_OVERLOAD``. Failed rows appear
here even when omitted from ``_delta.csv`` / ``_flat.csv``.

``<outputFile>_violations.csv`` (every ``outputFormat``) — streamed row
per branch/voltage violation. Columns: ``event_idx, contingency, type,
element, mva_or_vpu, rate_or_limit, loading_percent, base_mva, delta,
severity``. Emitted when ``loading_percent > violationSeverityThreshold
× 100`` (branch) or ``v_pu`` outside ``[minVoltage, maxVoltage]``
(voltage). ``loading_percent`` divides by the rate picked under
``contingencyRating``. ``severity = "critical"`` when branch loading
≥ 105 % or |Δv| ≥ 0.05 pu, else ``"warning"``.

``<outputFile>_summary.json`` (every ``outputFormat``) — end-of-run
aggregate: ``total_contingencies``, ``converged``, ``diverged =
total − converged``, plus the per-status split ``islanded``,
``no_slack``, ``solver_diverged``, ``slack_overload`` (sum to
``diverged``). Also ``with_branch_violation`` / ``with_voltage_violation``
/ ``violation_rows``; single-element extremes ``worst_loading`` /
``worst_voltage_low`` / ``worst_voltage_high`` (``null`` if none);
name arrays ``contingencies_with_branch_violation`` and
``contingencies_with_voltage_violation``; a length-``topN`` ranked
``top_severe_contingencies`` list (violated cases first by severity,
non-violated after by composite PI); and an echo of the run config
(``contingency_rating``, ``voltage_limit_low``/``high``,
``severity_threshold``).

When monitor filters are active the data-row count of ``_delta.csv`` /
``_flat.csv`` equals ``|monitored branches| × |converged
contingencies|``.

Output ordering
~~~~~~~~~~~~~~~

Rows are not sorted by contingency. The driver streams each MPI rank's
results to a per-rank ``.part`` file and rank 0 concatenates them in
rank order, so the final file is grouped by rank and ordered by
completion within each rank. Column 1 (``event_idx``) preserves
input-deck order — sort downstream if needed::

   ( head -1 my_run_delta.csv && \
     tail -n +2 my_run_delta.csv | sort -t, -k1,1n ) > my_run_delta.sorted.csv

Aggregated ``.txt`` outputs (``writeStats=true``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When ``writeStats=true`` the driver also writes the StatBlock summary
files: ``vmag.txt``, ``vmag_mm.txt``, ``vang.txt``, ``vang_mm.txt``,
``pgen.txt``, ``pgen_mm.txt``, ``qgen.txt``, ``qgen_mm.txt``,
``pflow.txt``, ``pflow_mm.txt``, ``qflow.txt``, ``qflow_mm.txt``,
``perf_mm.txt``, ``perf_sum.txt``, ``line_flt_cnt.txt``. The file
``pq_change_cnt.txt`` is added when ``qlim=true``. These hold per-element
statistics (mean / RMS / min / max) **across all contingencies**, not
per-contingency rows; their layout is described in the application
``README.md``. Set ``writeStats=false`` to skip them when CSV output is
sufficient.

Contingency XML aliases
~~~~~~~~~~~~~~~~~~~~~~~

The element names in a contingency XML accept PSS/E-aligned aliases:

+----------------------------+------------+-------------------------------+
| Element                    | Alias      | Holds                         |
+============================+============+===============================+
| ``contingencyLineNames``   | ``CKT``    | Branch circuit ID (PSS/E      |
|                            |            | ``CKT`` field).               |
+----------------------------+------------+-------------------------------+
| ``contingencyGenerators``  | ``GenID``  | Generator ID (PSS/E ``ID``    |
|                            |            | field).                       |
+----------------------------+------------+-------------------------------+

Either name parses; mix-and-match within a file is fine.

Advanced behavior
~~~~~~~~~~~~~~~~~

The driver also includes automatic slack-bus transfer when the slack
generator is the tripped element, a slack-capacity check after each
solve, and island / lone-bus detection when a branch trip splits the
network. The underlying mechanisms (``checkAndTransferSlack``,
``restoreSlack``, ``checkSlackCapacity``, ``getIslandCount``,
``hasLoneBus``) are exposed by ``PFAppModule`` and described above.

State Estimation Module
-----------------------

The state estimation module can be used to set up and run a state
estimation calculation. It does not have the extra functions that the
power flow module contains for supporting contingency analysis, so the
interface is a bit smaller. In addition to a standard network
configuration file, the state estimation calculation needs a second file
consisting of measurements. This file has the format

::

       <Measurements>
         <Measurement>
           <Type>VM</Type>
           <Bus>1</Bus>
           <Value>1.0600</Value>
           <Deviation>0.0050</Deviation>
         </Measurement>
         <Measurement>
           <Type>PIJ</Type>
           <FromBus>1</FromBus>
           <ToBus>2</ToBus>
           <CKT>BL</CKT>
           <Value>1.5688</Value>
           <Deviation>0.0100</Deviation>
         </Measurement>
         <Measurement>
           <Type>QIJ</Type>
           <FromBus>1</FromBus>
           <ToBus>2</ToBus>
           <CKT>BL</CKT>
           <Value>-0.2040</Value>
           <Deviation>0.0100</Deviation>
         </Measurement>
         <Measurement>
           <Type>PI</Type>
           <Bus>1</Bus>
           <Value>2.3240</Value>
           <Deviation>0.0100</Deviation>
         </Measurement>
         <Measurement>
           <Type>QI</Type>
           <Bus>1</Bus>
           <Value>-0.1690</Value>
           <Deviation>0.0100</Deviation>
         </Measurement>
         <Measurement>
           <Type>VA</Type>
           <Bus>1</Bus>
           <Value>0.0000</Value>
           <Deviation>0.0100</Deviation>
         </Measurement>
         <Measurement>
           <Type>IIJ</Type>
           <FromBus>1</FromBus>
           <ToBus>2</ToBus>
           <CKT>BL</CKT>
           <Value>1.5920</Value>
           <Deviation>0.0100</Deviation>
         </Measurement>
         <Measurement>
           <Type>IJI</Type>
           <FromBus>2</FromBus>
           <ToBus>1</ToBus>
           <CKT>BL</CKT>
           <Value>1.5915</Value>
           <Deviation>0.0100</Deviation>
         </Measurement>
       </Measurements>

The supported measurement types are ``VM`` (voltage magnitude),
``VA`` (voltage angle), ``PI`` (real power injection),
``QI`` (reactive power injection), ``PIJ`` (real power flow),
``QIJ`` (reactive power flow), ``IIJ`` (current magnitude,
from-bus direction), and ``IJI`` (current magnitude, to-bus
direction). Measurements can appear on any element of
the network and multiple measurements are allowed on each element. The
state estimation module does not have any error checking ability to
determine if there are sufficient measurements to guarantee solvability,
if not enough measurements are available then the calculation will
simply crash or fail to converge. An example ``State_estimation`` input block is shown below:

::

   <?xml version="1.0" encoding="utf-8"?>
   <Configuration>
     <State_estimation>
       <networkConfiguration> IEEE14.raw </networkConfiguration>
       <measurementList> IEEE14_meas.xml </measurementList>
       <maxIteration>20</maxIteration>
       <tolerance>1.0e-6</tolerance>
       <dampingFactor>1.0</dampingFactor>
       <maxBadDataIterations>5</maxBadDataIterations>
       <badDataThreshold>3.0</badDataThreshold>
       <diagnosticOutputLevel>standard</diagnosticOutputLevel>
       <useVoltageConstraints>false</useVoltageConstraints>
       <minVoltage>0.9</minVoltage>
       <maxVoltage>1.1</maxVoltage>
       <LinearSolver>
         <PETScOptions>
           -ksp_type richardson
           -pc_type lu
           -pc_factor_mat_solver_type superlu_dist
           -ksp_max_it 1
         </PETScOptions>
       </LinearSolver>
     </State_estimation>
   </Configuration>

The ``dampingFactor`` parameter (default 1.0) scales the Newton-Raphson
update step to improve convergence stability. The
``maxBadDataIterations`` and ``badDataThreshold`` parameters control
automatic bad data detection and elimination. The
``diagnosticOutputLevel`` can be set to ``basic``, ``standard``, or
``detailed`` to control the verbosity of solver output. The
``useVoltageConstraints`` flag enables virtual voltage limit
measurements that constrain bus voltages between ``minVoltage``
(default 0.9) and ``maxVoltage`` (default 1.1) per unit.

The state estimation module is
represented by the ``SEAppModule`` class which is in the
``gridpack::state_estimation`` namespace. The ``gridpack.hpp``
file contains a definition for the state estimation network
``SENetwork``. After instantiating an ``SEAppModule`` object and
a shared pointer to an ``SENetwork``, the state estimation
calculation can read in an external network configuration file using the
command

::

   void readNetwork(boost::shared_ptr<SENetwork> &network,
       gridpack::utility::Configuration *config)

The ``Configuration`` object should already be pointing at an open
file containing a ``State_estimation`` block. Inside the
``State_estimation`` block there should be a
``networkConfiguration`` field containing the name of the network
configuration file. The file name is parsed directly inside the
``readNetwork`` method and does not need to be handled by the user.
Alternatively, the ``SENetwork`` object may have already been cloned
from an existing network and therefore there is no need to read in the
configuration from an external file and partition it across processors.
In this case, the ``SEAppModule`` can be assigned the network using
the command

::

   void setNetwork(boost::shared_ptr<SENetwork> &network,
       gridpack::utility::Configuration *config)

This function just assigns the existing network to an internal pointer,
as well as a pointer to the input file. It is much more efficient than
reading in the network configuration file, if the network already
exists. This can occur when different types of calculations are being
chained together.

Once a network is in place and has been properly distributed, the
measurements can be read in by calling the function

::

   void readMeasurements()

The name of the measurement file is in the input deck and a pointer to
this file has already been internally cached in the ``SEAppModule``
when the network was assigned. The measurement file name is stored in
the ``measurementList`` field within the ``State_estimation``
block.

Alternatively, measurements can be provided programmatically using

::

   void setMeasurements(std::vector<Measurement> &measurements)

where each ``Measurement`` struct contains the type, bus or branch
indices, measured value, and standard deviation. This is useful when
measurements are generated by another module or read from a non-XML
source.

The network object can be initialized and the exchange buffers set up by
calling the

::

   void initialize()

method followed by

::

   void solve()

to obtain the solution to the system. The solver includes automatic bad
data detection. After each Newton-Raphson convergence, a chi-square test
is performed on the weighted residuals. If the chi-square test fails,
the measurement with the largest normalized residual is identified
and removed. The solver then re-converges with the remaining
measurements. This process repeats up to ``maxBadDataIterations``
times (default 5). Measurements with normalized residuals exceeding
``badDataThreshold`` (default 3.0) are flagged as bad data.

Results can be written out to standard out using the method

::

   void write()

This function will write out the voltage magnitude and phase angle for
each bus and the real and imaginary parts of the reactive power for each
branch. In addition, it will print out a comparison of the calculated
value and the original measured value for all measurements.

Finally, the results of the state estimation calculation can be saved to
the ``DataCollection`` object assigned to the buses by calling the

::

   void saveData()

method. The voltage magnitude and phase angle are stored as the
variables ``BUS_SE_VMAG`` and ``BUS_SE_VANG`` and the generator
parameters are stored as the indexed variables
``GENERATOR_SE_PGEN[i]`` and ``GENERATOR_SE_QGEN[i]``, where
``i`` runs over the set of generators on the bus.

The convergence status of the most recent solve can be queried using

::

   bool hasConverged() const

This returns ``true`` if the Newton-Raphson iteration converged within
the specified tolerance before reaching the maximum number of iterations.

Dynamic Simulation Module using Full Y-Matrix
---------------------------------------------

GridPACK supplies a dynamic simulation module that integrates the
equations of motion using an algorithm based on inversion of the full
Y-matrix. This module has been designed to enable the addition of
generator models that extend beyond the classical generator. It also
supports exciters, governors, power system stabilizers (PSS), relays,
dynamic loads, and inverter-based resource (IBR) / renewable energy models.
Models that are currently available include Generators:

::

   GENCLS
   GENSAL
   GENROU

IBR/Renewable Models:

::

   REGCA1
   REGCB1
   REGCC1
   GDFORM
   REGFMA1
   REECA1
   REPCA1

The ``EPRIA1`` model is also available when GridPACK is built with
``ENABLE_EPRI_IBR_MODEL``.

``REGFMA1`` is a droop-controlled grid-forming inverter model following
the WECC REGFM_A1 specification.

Exciters:

::

   EXDC1
   EXDC2
   ESDC2A
   ESST1A
   ESST3A
   ESST4B
   EXST1
   IEEET1
   IEEET2
   SCRX
   SEXS

Governors:

::

   WSIEG1
   WSHYGP
   TGOV1
   GAST
   HYGOV
   GGOV1
   IEEEG1
   IEESGO

Power System Stabilizers:

::

   PSSSIM
   IEEEST
   ST2CUT
   STAB2A

Wind Turbine Models:

::

   WTARA1
   WTDTA1
   WTPTA1
   WTTQA1

Relays:

::

   LVSHBL
   FRQTPAT
   DISTR1

Dynamic Loads:

::

   ACMTBLU1
   IEEL
   MOTORW
   CIM6BL

The full Y-matrix implementation of dynamic simulation is represented by
the ``DSFullApp`` class and the ``DSFullNetwork``, both of which
reside in the ``gridpack::dynamic_simulation`` namespace. The
dynamic simulation module uses an input deck of the form

::

   <?xml version="1.0" encoding="utf-8"?>
   <Configuration>
     <Dynamic_simulation>
       <networkConfiguration>IEEE_145.raw</networkConfiguration>
       <generatorParameters>IEEE_145.dyr</generatorParameters>
       <simulationTime>30</simulationTime>
       <timeStep>0.005</timeStep>
       <equilibriumInit>true</equilibriumInit>
       <Events>
         <faultEvent>
           <beginFault> 2.00</beginFault>
           <endFault>   2.05</endFault>
           <faultBranch>6 7</faultBranch>
           <timeStep>   0.005</timeStep>
         </faultEvent>
       </Events>
       <generatorWatch>
         <generator>
          <busID> 60 </busID>
          <generatorID> 1 </generatorID>
         </generator>
         <generator>
          <busID> 112 </busID>
          <generatorID> 1 </generatorID>
         </generator>
       </generatorWatch>
       <generatorWatchFrequency> 1 </generatorWatchFrequency>
       <generatorWatchFileName>gen_watch.csv</generatorWatchFileName>
       <LinearMatrixSolver>
         <PETScOptions>
           -ksp_atol 1.0e-18
           -ksp_rtol 1.0e-10
           -ksp_monitor
           -ksp_max_it 200
           -ksp_view
         </PETScOptions>
       </LinearMatrixSolver>
     </Dynamic_simulation>
   </Configuration>

The input for dynamic simulation module is contained in the
``Dynamic_simulation`` block. The optional ``equilibriumInit`` parameter
(default false) enables iterative Norton-equivalent initialization, which
eliminates pre-fault drift by rebalancing generator states (Pmech=Telec,
Efd=LadIfd) at the network-solved voltage. This is recommended for all
simulations.

Two features are important, the blocks
describing faults and the blocks describing monitored generators. Faults
are described in the ``Event``\ s block. The code currently only
handles faults on branches. Inside the ``Events`` block are
individual faults, described by a ``faultEvent`` block. Multiple
``faultEvent`` blocks can be contained within the
``Events`` block. As will be described below, it is possible
for the faults to be listed in a separate file. This can be convenient
for describing a task-based calculation that may contain a lot of
faults. The parameters describing the fault include the time (in
seconds) that the fault is initiated, the time that it is terminated,
the timestep used while integrating the fault and the indices of the two
buses at either end of the fault branch.

When running a dynamic simulation, it is generally desirable to monitor
the behavior of a few generators in the system and this can be done by
setting generator watch parameters. The ``generatorWatch`` block
specifies which generators are to be monitored. Each generator is
described within a ``generator`` block that contains the index of
the bus that the generator is located on and the character string ID of
the generator. The results of monitoring the generator are written to
the file listed in the ``generatorWatchFileName`` field and the
frequency for storing generator parameters in this file is set in the
``generatorWatchFrequency`` field. This parameter describes the time
step interval for writing results (an integer), not the actual time
interval.

Before using the dynamic simulation module, a network needs to be
instantiated outside the ``DSFullApp`` and then passed into the
module. If the module itself is going to read and partition a network,
then it should use the function

::

   void readNetwork(boost::shared_ptr<DSFullNetwork> &network,
     gridpack::utility::Configuration *config,
     const char *otherfile = NULL)

The ``Configuration`` object should already be pointing to an input
deck with a ``Dynamic_simulation`` block that specifies the network
configuration file. The optional ``otherfile`` argument in
``readNetwork`` can be used to overwrite the
``networkConfiguration`` field in the input deck with a different
filename. This capability has proven useful in some contingency
applications where multiple PSS/E files need to be read.

Alternatively, a distributed network may already exist (it may have been
cloned from another calculation). In that case, the function

::

   void setNetwork(boost::shared_ptr<DSFullNetwork> &network,
     gridpack::utility::Configuration *config)

can be used to assign an internal pointer to the network. Again, the
``Configuration`` object should already be pointing to an input
file.

Additional generator parameters can be assigned to the generators by
calling the function

::

   void readGenerators()

This function opens the file specified in the
``generatorParameters`` field in the input file and reads the
additional generator parameters. The file is assumed to correspond to
the PSS/E .dyr format. The devices listed at the start of this section
can be included in this file.

After setting up the network and reading in generator parameters, the
module can be initialized by calling

::

   void initialize()

This sets up internal parameters and initializes the network so that it
is ready for calculations.

A list of faults can be generated from the input file by calling

::

   std::vector<gridpack::dynamic_simulation::DSFullBranch::Event>
     getEvents(gridpack::utility::Configuration::CursorPtr cursor)

If the cursor variable is pointed at a ``Dynamic_simulation`` block
inside the input file (as in the example input block above) then this
function will return a list of faults from the input deck. However, it
is also possible that the cursor could be pointed to the contents of
another file. As long as it is pointed to a block containing a
``Events`` block, this function will return a list of faults. This
allows users to declare a large list of faults in a separate file and
then access the list by including the external file name as a parameter
in the input deck of their application.

The monitoring of generators specified in the input deck can be set up
by calling

::

   void setGeneratorWatch()

This will guarantee that all generators specified in the input deck are
monitored and that the results are written out to the specified file. If
this function is not called, the generator watch parameters in the input
file are ignored. Simulations can be run using the function

::

   void solve(gridpack::dynamic_simulation::DSFullBranch::Event fault)

Some additional results can be written at the end of the simulation
using the function

::

   void write(const char *signal)

The signal parameter can be used to control which results are written
out. This function currently does not support any output. All output
results are controlled using the generator watch parameters.

Some additional functions can be used to control where output generated
during the course of a simulation is directed. The following two
functions can be used to direct output from the ``write`` function
to a file

::

   void open(const char* filename)

   void close()

The function

::

   void print(const char* buf)

can be used to print out a string to standard out. If the ``open``
function has been used to open a file, then the output is directed to
the file. This function is equivalent to the ``header`` convenience
function in the serial IO classes. Additional functions can be used to
stored data from the generator watch variables. These can be used to
save the time series data from a simulation in a collection of vectors.
The application can then use these series in whatever way it wants.
There are four functions that enable this capility. The first is

::

   void saveTimeSeries(bool flag)

This function must be called with the argument set to “true” in order
for the time series data to be saved. Otherwise it is only written to
output and no data is saved between time steps. The second function can
be called after the solve function has been called and the simulation is
completed. It returns a vector of time series

::

   std::vector<strd::vector<double> > getGeneratorTimeSeries()

This function returns a vector containing the time series data for all
the watched generators located *on this processor* (generators on buses
owned by neighboring processors are not included).

To find out which variables are actually in the list returned by
``getGeneratorTimeSeries`` requires the remaining two functions. The
function

::

   void getListWatchedGenerators(std::vector<int> &bus_ids,
       std::vector<std::string> &gen_ids)

returns a list of the bus IDs and 2-character generator tags for all
monitored generators. In particular, it assigns and ordering to these
generators that is used by function

::

   std::vector<int> getTimeSeriesMap()

This function returns a map between the elements in the list of time
series returned by ``getGeneratorTimeSeries`` and the generators
that those time series correspond to. For example suppose the time
series list has four elements in it that happen to correspond to two
generators on processor. There are a total of six monitored generators
in the system. The vectors returned by ``getListWatchedGenerators``
have length six, the vector returned by ``getTimeSeriesMap`` has
length four. The value in the map vector for the corresponding element
in the time series vector points to the location of the bus index and
generator tag for that time series variable in the lists returned by
``getListWatchedGenerators``. This still leaves it up to the user to
identify the actual variable being watched within the generator. In this
example there are four variables that are watched but only two
generators. Currently, the generator watch capability only watches the
rotor speed and rotor angle of each generator. The first time series is
the speed and the second time series is the angle.

Recently, a new function was added to the dynamic simulation application
module that allows users to update the data collection objects for all
branches an buses in the system with current values from the simulation.
A call to the function

::

   void updateData()

will take the current values of internal variables and copy them to
values in each network component’s data collection object. These values
can then be accessed by querying the appropriate data collection object.
Variables in the dictionary header files that are updated by the
``updateData`` function have ``_CURRENT`` appended to their
names. At present, the following variables are updated by
``updateData``

::

   BRANCH_FROM_P_CURRENT
   BRANCH_FROM_Q_CURRENT
   BRANCH_TO_P_CURRENT
   BRANCH_TO_Q_CURRENT
   BRANCH_IRFLOW_CURRENT
   BRANCH_IIFLOW_CURRENT
   BUS_VMAG_CURRENT
   GENERATOR_PG_CURRENT
   GENERATOR_QG_CURRENT
   GENERATOR_VS_CURRENT
   GENERATOR_NOM_POWER_CURRENT
   LOAD_PL_CURRENT
   LOAD_QL_CURRENT

Other values can be added as needed. Each of the network components, bus
and branch, have an ``updateData`` method that is called by the
``updateData`` method in the ``DSFullApp`` class. This also
extends to the base model classes, which all have default
implementations that do not do anything. Individual models can overwrite
the default implementations to add data to the data collection objects
when the update function is called.

Kalman Filter
-------------

GridPACK includes a Kalman filter module that can be used for dynamic
state estimation. The Kalman filter relies heavily on parallel matrix
multiplies that are not currently very high performing, so users will
probably find this module too slow for large grids. However, we include
it for users interested in exploring the use of Kalman filters in
smaller applications. We hope to improve performance in future releases.

The current implementation of the Kalman filter only supports classical
generators. These are described in a PSS/E .dyr formatted file. The
network itself can be described using a standard PSS/E .raw file. In
addition to the .raw and .dyr files, users need to supply times series
data for the voltage magnitude and voltage phase angle on all buses.
These are stored as .csv files. The format for both the voltage
magnitude and phase angle files is

::

   t-3001,  Bus-1,  Bus-2,{\dots}
   0.0,    -0.001, -0.135,{\dots}
   0.1,    -0.001, -0.135,{\dots}
   0.2,    -0.001, -0.135,{\dots}
    :

All entries on the same lines are separated by commas. The first row
contains the name of all columns. The first column is time and has a
name of the form ``t-xxx,`` where ``xxx`` is an integer
representing the number of time steps in the file. The number of rows in
the file corresponds to ``xxx+1`` (the extra row is the first line
with the column names). The number of columns is equal to the number of
buses in the file plus one (the extra column contains the times). After
the first column, the remaining names are all of the form
``Bus-xxx``, where ``xxx`` is an integer representing the bus
ID. The remaining rows contain the time of the measurement and the value
for the measurement on each of the buses.

The input file for the Kalman filter module used both for a dynamic
simulation as well as input that is unique to the Kalman filter module.
The dynamic simulation parameters that are used include

::

     <Dynamic_simulation>
       <simulationTime>3</simulationTime>
       <timeStep>0.01</timeStep>
       <!-- = 1 Fault Event is known; 
            = 0 Fault event is unknown, switch is skipped. 
       -->
       <KnownFault> 1 </KnownFault>
       <TimeOffset> 0 </TimeOffset> <!--skip initial measurement data -->
       <Events>
         <faultEvent>
           <beginFault> 1 </beginFault>
           <endFault>   1.1</endFault>
           <faultBranch>6 7</faultBranch>
           <timeStep>   0.01</timeStep>
         </faultEvent>
       </Events>
     </Dynamic_simulation>

The fault used in the simulation is specified using the same
``Events`` block as for dynamic simulation. If the Kalman filter
simulation is not being initialized from another calculation, the
``networkConfiguration`` field can also be added. The KnowFault and
TimeOffset parameters are unique to the Kalman filter application and
control whether the fault is considered to be a know event and whether
all the time series data should be used in the analysis.

The Kalman filter block consists of the fields

::

     <Kalman_filter<
       <KalmanAngData>IEEE14_Kalman_input_ang.csv</KalmanAngData>
       <KalmanMagData>IEEE14_Kalman_input_mag.csv</KalmanMagData>
       <generatorParameters>IEEE14_classicGen.dyr</generatorParameters>
       <ensembleSize>21</ensembleSize>
       <gaussianWidth>1e-2</gaussianWidth>
       <noiseScale>1e-4</noiseScale>
       <randomSeed>931316785</randomSeed>
       <maxSteps>3000</maxSteps>
       <LinearSolver>
         <PETScOptions>
           -ksp_view
           -ksp_type richardson
           -pc_type lu
           -pc_factor_mat_solver_package superlu_dist
           -ksp_max_it 1
         </PETScOptions>
       </LinearSolver>
     </Kalman_filter>

The ``KalmanAngData`` and ``KalmanMagData`` fields specify the
locations of the files containing the time series data for the voltage
magnitude and phase angle. The .dyr file containing the generator
parameters (classical generators only) is specified in the
``generatorParameters`` field. Additional Kalman filter parameters
include

#. ``ensembleSize``: The number of random ensembles generated for
   the Kalman filter calculation.

#. ``gaussianWidth``:

#. ``noiseScale``:

#. ``randomSeed``: This is an arbitrary integer used to seed the
   GridPACK random number generator.

#. ``maxSteps``: this parameter can be used to control the number of
   steps simulated. If the number of steps is smaller than the number of
   steps in the time series data files, then only the number of steps
   set by ``maxSteps`` will be simulated.

The Kalman filter also needs to make use of linear solvers and the type
of solver and its parameters can be specified in this block as well.

The Kalman filter module is represented by the ``KalmanApp`` class
and the ``KalmanNetwork``, both of which are in the
``gridpack::kalman_filter`` namespace. At present there are only a
few functions in this class, more will probably be added as we develop
this module further. Apart from the constructor and destructor, the
``KalmanApp`` class has a method for reading in a network from a
PSS/E formatted file and partitioning it among processors

::

   void readNetwork(boost::shared_ptr<KalmanNetwork> &network,
       gridpack::utility::Configuration *config)

If the network already exists, then it can be applied to an existing
KalmanApp object using the function

::

   void readNetwork(boost::shared_ptr<KalmanNetwork> &network,
       gridpack::utility::Configuration *config)

The application can be initialized by calling the function

::

   void initialize()

This function will read in the files containing the time series data for
the voltage magnitude and phase angles and will set update configure the
calculation based on the parameters in the input file. The simulation is
run and output generated using

::

   void solve()

The values of the rotor speed and rotor angle for all generators will be
written to the files ``omega.dat`` and ``delta.dat`` after this
simulation is run.

Electromagnetic Transient (EMT) Module
---------------------------------------

GridPACK includes an electromagnetic transient (EMT) simulation module
that provides detailed modeling of power system components using 
differential-algebraic equations (DAE). The EMT module is designed for
high-fidelity analysis of electromagnetic phenomena in power systems,
making it suitable for studying transient behavior, fault analysis,
and detailed generator dynamics. The module uses a DAE solver to integrate
the system equations and supports both explicit and implicit integration
algorithms for machine models.

The EMT module is represented by the ``Emt`` class and uses an 
``EmtNetwork`` for the network representation. Unlike other GridPACK
modules, the EMT class does not use a specific namespace but is designed
to work seamlessly with the GridPACK framework. The module requires both
a power flow initialization phase and an EMT simulation phase.

The EMT module uses an input deck with both ``Powerflow`` and ``EMT`` 
configuration blocks. An example input file has the form:

::

   <?xml version="1.0" encoding="utf-8"?>
   <Configuration>
     <Powerflow>
       <networkConfiguration_v34>case9_PSCAD.raw</networkConfiguration_v34>
       <maxIteration>50</maxIteration>
       <tolerance>1.0e-6</tolerance>
       <LinearSolver>
         <PETScOptions>
           -ksp_type richardson
           -pc_type lu
           -pc_factor_mat_solver_package superlu_dist
           -ksp_max_it 1
         </PETScOptions>
       </LinearSolver>
       <UseNonLinear>false</UseNonLinear>
     </Powerflow>
     <EMT>
       <generatorParameters>case9_GAST.dyr</generatorParameters>
       <machineIntegrationType>EXPLICIT</machineIntegrationType>
       <simulationTime>2.0</simulationTime>
       <timeStep>0.00005</timeStep>
       <Events>
         <BusFault>
           <begin>0.5</begin>
           <end>0.6</end>
           <bus>9</bus>
           <type>ThreePhase</type>
           <Ron>0.011</Ron>
           <Rgnd>0.0</Rgnd>
         </BusFault>
       </Events>
       <Monitors>
         <Generator>
           <bus>1</bus>
           <id>1</id>
         </Generator>
         <Generator>
           <bus>2</bus>
           <id>1</id>
         </Generator>
       </Monitors>
       <DAESolver>
         <PETScPrefix>emt_</PETScPrefix>
       </DAESolver>
     </EMT>
   </Configuration>

The ``Powerflow`` block contains standard power flow configuration 
parameters including the network configuration file and solver options.
The ``EMT`` block contains parameters specific to the EMT simulation:

#. ``generatorParameters``: PSS/E .dyr format file containing additional
   generator model parameters for the EMT simulation.

#. ``machineIntegrationType``: Integration algorithm for machine models.
   Can be set to ``EXPLICIT`` or ``IMPLICIT`` (default).

#. ``simulationTime``: Total simulation time in seconds.

#. ``timeStep``: Integration time step for the EMT simulation.

#. ``Events``: Block containing fault events and other disturbances.
   Currently supports ``BusFault`` events with parameters for fault
   timing, location, type, and impedance values.

#. ``Monitors``: Specifies generators to be monitored during simulation.
   Output from monitored generators is written to files.

#. ``DAESolver``: Configuration parameters for the DAE solver.

After instantiating an ``Emt`` object with a communicator, the EMT
calculation can be set up using the following sequence of function calls.
First, the configuration file must be assigned:

::

   void setconfigurationfile(const char* configfile)

This function assigns the input configuration file to the EMT object.
Next, an initial power flow calculation must be performed to establish
the system operating point:

::

   void solvepowerflow()

This function reads the network configuration, performs the power flow
calculation, and initializes the system state variables. The power flow
solution provides the initial conditions for the EMT simulation.

After the power flow is complete, the EMT-specific setup can be performed:

::

   void setup()

This function reads the generator parameters from the .dyr file, sets up
the EMT network components, initializes the DAE solver, and prepares the
system for transient simulation. The setup process includes creating
vector and matrix mappers for the DAE system and configuring event
handlers for faults and other disturbances.

The EMT simulation can then be executed using:

::

   void solve()

This function integrates the DAE system over the specified simulation time,
handling events such as faults and writing output for monitored generators.
The integration uses the configured time step and machine integration type.

The EMT module supports several types of events that can be specified in
the input configuration:

#. **Bus Faults**: Three-phase faults applied to specified buses with
   configurable fault resistance and grounding resistance.

#. **Generator Monitoring**: Real-time monitoring and output of generator
   state variables during simulation.

The module uses a DAE solver framework that supports various integration
algorithms and can handle stiff differential equations typical in power
system transient analysis. The solver configuration can be controlled
through PETSc options in the ``DAESolver`` block.

Output from the EMT simulation includes time series data for monitored
generators, which is written to files during the simulation. The output
frequency and monitored variables can be controlled through the input
configuration.

The EMT module is designed to work in conjunction with other GridPACK
modules, particularly power flow, to provide comprehensive power system
analysis capabilities. The module can be used as a starting point for
developing specialized EMT applications or can be integrated into larger
multi-physics simulations.
