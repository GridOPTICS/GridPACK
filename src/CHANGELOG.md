# Change Log
The format is based on [Keep a Changelog](http://keepachangelog.com/).

All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).

This project follows the [Gitflow Workflow
model](https://www.atlassian.com/git/tutorials/comparing-workflows/gitflow-workflow).

## [Unreleased]
The Unreleased section will be empty for tagged releases. Unreleased
functionality appears in the develop branch.

## [3.7.0]
- Added
  - Power Flow
    - Switched shunt control (discrete MODSW=1 and continuous MODSW=2 modes)
      with per-block N/B stepping, cycle detection, and lockout
    - LTC (load tap changer) transformer control with tap direction
      dependent on controlled bus side (from-bus vs to-bus)
    - Area interchange control for power flow
    - Auto-detect PSS/E version from RAW file header
    - IREG remote voltage regulation via PV bus swap: for generators with
      IREG != 0, the remote bus becomes PV (V = VS) and the generator bus
      becomes PQ, matching PSS/E and PowerWorld behavior
    - Q-limit deadband (qlimDeadband XML parameter) to prevent PV/PQ
      cycling at marginal Q violations
    - NR step damping for improved convergence robustness
    - Test cases for switched shunt, LTC, and area interchange controls
  - Contingency Analysis
    - Monitor filters: monitorBranchesFile allowlist, monitorAreas, and
      monitorKvMin/Max, applied to every output format
    - Output formats csv, csv_flat, csv_delta and json, plus per-run
      _violations.csv, _convergence.csv and _summary.json with performance
      index ranking and worst-contingency rosters
    - contingencyRating (A|B|C, default A) as the loading% denominator
      across all outputs
  - Dynamic Simulation
    - REGFMA1 grid-forming model; QVflag initialization option
    - Exciter models EXDC2, ESDC2A, ESST3A, EXST1, IEEET2 and SCRX;
      governor IEESGO; PSS models IEEEST, ST2CUT and STAB2A
  - EMT
    - Ported dynamic simulation model fixes (saturation, auto-corrections,
      limit guards) to the EMT layer; updated SEXS and REECA1 ports
  - Python
    - updateData stores voltage phase angle
  - Documentation
    - User manual updated with switched shunt, LTC and contingency analysis
      output/filter sections
- Fixed
  - Power Flow
    - Star bus output filtering: synthetic buses from 3-winding transformers
      are excluded from power flow output
    - SHUNT_SWREM/SWREG name mismatch: v34 parser stored remote regulation
      bus under SWREG but PF component read SWREM; added fallback
    - LTC tap direction: reversed tap adjustment for to-bus controlled
      transformers (COD1=1, CONT1=to-bus)
    - 3-winding transformer ratings in v34/v35 parsers: winding 2 and 3
      ratings were written to data1 instead of data2/data3
    - 3-winding transformer star bus voltage initialization from
      VMSTAR/ANSTAR fields
    - Warm-start Q-limit clamp picked the wrong limit for generators whose
      scheduled QG disagreed with the saved voltage
    - Phantom flow reported on out-of-service branches
  - Contingency Analysis
    - StatBlock crash ("wrong dimension specified") when a monitor filter
      left no buses, generators or branches
  - Parser
    - 2-winding transformer ratings A/B/C not populated from RATE1-3 in
      v34/v35
    - PSS/E v34 switched shunt parser: PTI34 incorrectly used v35 parser
      (SwitchedShuntParser35) which assumed ID and NREG fields not present
      in v34 format, causing BINIT to read B1 instead (e.g., 100 vs 400
      Mvar), producing large voltage errors at switched shunt buses
  - Data
    - IEEE118 RAW files: base kV corrected on 56 buses (54 were 1 kV,
      buses 57/58 mistyped); 106 buses at 138 kV, 12 at 345 kV
  - Python
    - Removed distutils import that broke on Python 3.12+
  - Build
    - Install missing parser_classes headers used by base_pti_parser

## [3.6]
- Added
  - Docker Support
    - Multi-architecture Docker builds (AMD/ARM) with CI/CD pipeline
    - Complete containerized GridPACK installation (OpenMPI, Boost, GA, PETSc)
    - Documentation updated to include Docker usage
  - State Estimation Enhancements
    - Jacobian optimization and performance profiling for large-scale systems
    - Bad data detection with chi-square test and configurable
      diagnosticOutputLevel parameter (basic/standard/detailed)
    - Sparse matrix support for improved computational efficiency
    - Added IIJ/IJI measurement types and enhanced voltage constraint options
    - HPC capability for state estimation
    - Chi-square convergence criteria in the outer loop
    - Comprehensive output reporting with convergence warnings
    - XML-configurable Newton-Raphson damping factor (dampingFactor, default 1.0)
    - Updated IEEE 14-bus and IEEE 118-bus test inputs
  - Python interface for state estimation including setMeasurements() API
    and Sphinx documentation for the Python state estimation interface.
  - Power Flow Improvements
    - Q-limit (QLIM) handling with PV-to-PQ bus switching and island detection
    - RMPCT-based reactive power distribution for multi-generator buses
    - Fixed shunt status check, Qg distribution when QLIM is reached, and
      flat/warm start 
    - Switched shunt support
    - ZIP load model with voltage-dependent load representation
    - Added Pinj, Qinj to power flow screen outputs
    - Per-iteration max mismatch reporting in Newton-Raphson solve
    - JSON and CSV export of power flow results via ResultsExporter
  - Contingency Analysis Enhancements
    - N-1 auto-generation feature for contingency analysis
    - Automatic slack bus transfer and capacity check
    - QLIM configuration support in CA input XML files (default: false)
    - PV->PQ warning messages in contingency output files
    - Generator ID display in contingency results
    - JSON and CSV export of contingency analysis results
  - Dynamic Simulation Enhancements
    - GENROU PSS integration
    - Iterative equilibrium initialization (XML: equilibriumInit, default false)
    - Support dynamic models (GENROU, GENSAL, Classical) for generators with
      negative PG.
  - Added build trigger when releasing a new version of GridPACK
- Changed
  - User manual updated for v3.6 including power flow, state estimation,
    dynamic simulation, and Python interface documentation
  - Migrated QLIM parameter to bool type in power flow and CA for consistency
  - Optimized CMake builds
- Fixed
  - Generator Models
    - GENROU: scaled quadratic saturation, iterative q-axis saturation init,
      LadIfd Se*Psidpp term, (1+omega) speed factor on internal voltage
    - GENSAL: unified saturation to PSS/E scaled quadratic, (1+omega) speed
      factor, predictor speed deviation, saturation initialization
    - PowerAngle output wrapping to [-180,180] in all generator models
    - Null pointer crash on unrecognized generator model in dsf_components
  - Exciter Models
    - ESST1A: regulator limits Vrmin/Vrmax to Vamin/Vamax, lead-lag bypass
    - ESST4B: Kim parser key from EXCITER_KPM to EXCITER_KIM
    - EXDC1: lead-lag using TA instead of TC
    - SEXS: EMAX default from 0.0 to 999.0
  - Governor Models
    - GGOV1: re-enabled parameter parsing, PIDIn undefined behavior, State 9
      ODE, lead-lag bypass/coefficients, corrector rate limit variable
    - WSIEG1: parser T3 field, NGV yin[] index, NGV/Db1/Db2 block types and
      conditional bypass for zero data
    - WSHYGP: NGV yin[] indices, gate servo feedback, abs/fabs
    - HYGOV: fallback setting R instead of r
    - PSSSIM: washout output divided by Tw
  - Renewable Energy Models
    - REECA1: VDL1/VDL2 swap, Thld state machine, timer formula, TIQ parser
    - REPCA1: FEMIN parser column, over-frequency droop sign handling
    - REGCA1: Iq_olim sign, P/Q swap in output
    - REGCB1/REGCC1: uninitialized Preal variable
  - Load Models
    - MOTORW: comma operator (C0=1,0 -> 1.0), hardcoded llr1 removed, slip
      interpolation
    - IEEL: normalization destroying valid coefficients
  - Framework / Infrastructure
    - ForceSerial off-by-one in setElementRange half-open range
    - Cblock crash when model time constant is zero
    - Contingency analysis multi-process deadlock in GA operations
    - 3-winding transformer parsing in PSS/E parsers
  - Power Flow
    - Warm start with qlim=false scenario
    - Incomplete printing of parallel lines in output
    - Trailing space in contingency output filenames
    - Qmin/Qmax storage

## [3.5]
- Added
  - A new application for electromagnetic transient simulation has been added to GridPACK. Its features are as follows:
    - Has typical models of synchronous generators (GENROU) and controls (exciter and turbine governor)
    - Includes Grid-Following (REGCA1) and Grid-Forming (REGCFM) inverter models
    - Transmission lines are modeled as lumped parameter lines
    - Loads are modeled as constant impedance loads
    - Uses .raw and .dyr as inputs.
    - Steady-state initialization
    - Fixed-step and variable time-stepping support
    - Parallelization through network partitioning
    - Number of test cases including IEEE-9bus, IEEE-39bus, Kundur 2-area network, and WECC-240 bus.
  - Added install_gridpack_deps.sh and install_gridpack.sh scripts to build
    libraries used by GridPACK and to build and install GridPACK based on those
    libraries. These scripts work most of the time, but may need to be
    customized for individual platforms.
  - GridPACK can be initialized from a user-supplied communicator (if using GA
    5.9 or greater) instead of forcing GridPACK to initialize on MPI_COMM_WORLD.
    This can be useful if integrating GridPACK with other applications that
    may be spawning off GridPACK simulations as individual tasks.
  - Added new version of dynamic simulation based on variable time-stepping
    algorithm. This can potentially run much faster than fixed timestep
    implementations since large timesteps can be used when the system is
    relatively stable or changes are occuring on long time scales.
  - Added a parser for Mat Power files to GridPACK. This only supports parsing
    of .baseMVA, .bus, .gen, .branch, .areas and .gencost blocks.
  - Added a singleton NoPrint object that can be used to suppress all external
    printing.
  - Added export modules for PSS/E v23, v33, and v34 formatted files. Not all
    data blocks and variables are supported so in general a file that is read
    in and then exported will have less data than the original file. In general,
    if a variable or block is not used by the current GridPACK applications,
    it is likely that it will not appear in the exported file.
  - Expanded the dynamic simulation application by added many new models. These
    include
    - Added a grid-forming inverter-based resource (IBR) model (e.g. gdform) and
      grid-following IBR models (e.g. regca1, regcb1, regcc1, epria1) and control
      model (reeca1).
    - Added block diagram components (e.g. DelayBlock, DelayBlockwithLimit,
      PIBlockwithLimit) and updated old generator, exciter and governor models by
      block diagram implementation.
    - Added several new exciter and governor models (e.g. gast, hygov, ieeet1, sexs,
      tgov1).
    - Implemented wind turbine control blocks (e.g. wtara1, wtdta1, wtpta1, wttqa1).
  - Added HADREC (Hierarchical Adaptive Dynamic Resilience Coordinator) module
    providing real-time grid monitoring and control capabilities including load
    shedding, line/generator tripping, wide-area control signals for PSS, and
    zone-based load/generation tracking with observation collection interfaces.
  - Completed Python interface to power flow, HADREC, dynamic simulation
    application modules.
  - Added electromagnetic transient (EMT) simulation capability supporting
    three-phase abc domain analysis with DAE-based time-domain simulations
    and abc-to-dq0 reference frame transformations.
  - Implemented variable timestepping algorithms in dynamic simulation and EMT
    modules, enabling adaptive timestep control for improved computational
    efficiency and numerical stability during system transients.
  - Updated parsers and dictionary to handle wind machines.
- Changed
  - The user manual has been moved to a Github ReadTheDocs location and is now
    available on the web. The previous PDF files are no longer supported.
  - The user manual now includes Python interface documentation.
  - The webpages at www.gridpack.pnl.gov based on Wikimedia are no longer
    supported. These webpages have been converted to markdown syntax and are
    now part of the GitHub repository. Documentation can be found by scrolling
    to the README.md section of the GitHub repository for GridPACK and following
    the links from there.
  - Updated Math Module DAE solver to support both real and complex datatypes
- Fixed
  - Fixed bug in PSS/e parsers so that names containing '\' character are not
    confused with comments.
  - Modified PSS/E parsers so that they can read files that use the convention
    "value1,value2,,,,value6,value7....". The missing values are assumed to be
    0. True PSS/E parsers may set these to a default value. If this is not 0,
    then these values will fail in GridPACK.
- Known Limitations
  - 3-winding transformers are not fully supported PSS/E version 33-36 parsers.
