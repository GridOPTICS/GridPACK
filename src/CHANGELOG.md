# Change Log
The format is based on [Keep a Changelog](http://keepachangelog.com/).

All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).

This project follows the [Gitflow Workflow
model](https://www.atlassian.com/git/tutorials/comparing-workflows/gitflow-workflow).

## [Unreleased]
The Unreleased section will be empty for tagged releases. Unreleased
functionality appears in the develop branch.

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
      flat/warm start (VM/VA)
    - Switched shunt support
    - ZIP load model with voltage-dependent load representation
    - Added Pinj, Qinj to power flow screen outputs
  - Contingency Analysis Enhancements
    - N-1 auto-generation feature for contingency analysis
    - Automatic slack bus transfer and capacity check
    - QLIM configuration support in CA input XML files (default: false)
    - PV->PQ warning messages in contingency output files
    - Generator ID display in contingency results
  - Added build trigger when releasing a new version of GridPACK
- Changed
  - User manual updated for v3.6 including power flow, state estimation,
    dynamic simulation, and Python interface documentation
  - Migrated QLIM parameter to bool type in power flow and CA for consistency
  - Optimized CMake builds
- Fixed
  - Fixed 3-winding transformer parsing in PSS/E parsers
  - Fixed warm start with qlim=false scenario
  - Fixed incomplete printing of many parallel lines in power flow output
  - Fixed trailing space in contingency output filenames
  - Fixed Qmin/Qmax storage

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
