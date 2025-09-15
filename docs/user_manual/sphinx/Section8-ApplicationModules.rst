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
file containing the network information. The network configuration file
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

   bool unsetContingency(Contingency &event)

can be used to set or unset a contingency. The call
``unsetContingency`` should only be called after calling
``setContingency`` and it should use the same ``event``
argument. After calling the ``unsetContingency`` method, the network
should have the same configuration as before calling the
``setContingency`` method.

The remaining calls in ``PFAppModule`` can be used to determine the
status of a network after solving a configuration with a contingency.
The functions

::

   bool checkVoltageViolations(double Vmin, double Vmax)

   bool checkVoltageViolations(int area, double Vmin, double Vmax)

can be used to check for a voltage violation anywhere in the system
where ``Vmin`` and ``Vmax`` are the minimum and maximum
allowable voltage excursions (per unit). The second function only checks
for violations on buses with the specified value of ``area``. These
functions are true if there are no voltage violations and return false
if a violation is found on one or more buses. It frequently turns out
that many networks have voltage violations even in the absence of any
contingencies and it is often desirable to ignore these violations. This
can be accomplished using the function

::

   void ignoreVoltageViolations(double Vmin, double Vmax)

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

The limits on the line are contained in parameters read in from the
network configuration file so these functions have no arguments
describing the line limits. The second function will only check for
violations on lines with the specified value of ``area``. Like
voltage violations, branches that display line overload violations that
are present even without contingencies can be ignored in the checks by
calling the function

::

   void ignoreLineOverloadViolations()

after running a calculation on the system without contingencies. These
settings can be cleared using the function

::

   void clearLineOverloadViolations()

Finally, the internal voltage variables that are used as the solution
variables in the power flow calculation can be reset to their original
values (specified in the network configuration file) by calling the
function

::

   void resetVoltages()

Again, this may be useful in contingency calculations where multiple
calculations are run on the same network and it is desirable that they
all start with the same initial condition.

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
           <Values>-0.2040</Value>
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
       </Measurements>

for the five types of measurements ``VM``, ``PIJ``, ``QIJ``,
``PI``, and ``PJ``. Measurements can appear on any element of
the network and multiple measurements are allowed on each element. The
state estimation module does not have any error checking ability to
determine if there are sufficient measurements to guarantee solvability,
if not enough measurements are available then the calculation will
simply crash or fail to converge. The state estimation module is
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

The network object can be initialized and the exchange buffers set up by
calling the

::

   void initialize()

method followed by

::

   void solve()

to obtain the solution to the system. Results can be written out to
standard out using the method

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

Dynamic Simulation Module using Full Y-Matrix
---------------------------------------------

GridPACK supplies a dynamic simulation module that integrates the
equations of motion using an algorithm based on inversion of the full
Y-matrix. This module has been designed to enable the addition of
generator models that extend beyond the classical generator. It also
supports exciters, governors, relays and dynamic loads. Models that are
currently available include Generators:

::

   GENCLS
   GENSAL
   GENROU

Exciters:

::

   EXDC1
   ESST1A

Governors:

::

   WSIEG1
   WSHYGP

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
``Dynamic_simulation`` block. Two features are important, the blocks
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
