
Power Flow
**********

The Python Power Flow interface mirrors closely the ``PFAppModule``
C++ interface.  See the `that module's documentation
<https://gridpack.readthedocs.io/en/latest/Section8-ApplicationModules.html#power-flow>`_
for most details.


Usage
=====

A minimal example:

.. code:: python

   import sys, os
   import gridpack
   import gridpack.powerflow

   env = gridpack.Environment()
   comm = gridpack.Communicator()
   config = gridpack.Configuration()
   config.open("input.xml", comm)

   pfapp = gridpack.powerflow.Powerflow()
   pfapp.readNetwork(config, -1)
   pfapp.initialize()
   pfapp.solve();
   pfapp.write();

   # try to force deletion order to avoid problems
   del pfapp
   del env


Reference
=========

GridPACK power flow module

**class gridpack.powerflow.Powerflow**

   Power flow application module

   **checkLineOverloadViolations(*args, **kwargs)**

      Overloaded function.

      1. checkLineOverloadViolations(self:
         gridpack.powerflow.Powerflow) -> bool

      Check to see if there are any line overload violations in the
      network

      1. checkLineOverloadViolations(self:
         gridpack.powerflow.Powerflow, arg0: int) -> bool

      Check to see if there are any line overload violations in an
      area

      1. checkLineOverloadViolations(self:
         gridpack.powerflow.Powerflow, arg0: list[int], arg1:
         list[int], arg2: list[str]) -> object

      Check for overload violations on specific lines

      Parameters:
            bus1 : list of int
               original index of “from” bus for branch

            bus2 : list of int
               original index of “to” bus for branch

            tags : list of int
               line IDs for individual lines

         violations true if violation detected on branch, false
         otherwise

      Returns:
         True if no violations found, list of flags for each line
         otherwise

   **checkQlimViolations(*args, **kwargs)**

      Overloaded function.

      1. checkQlimViolations(self: gridpack.powerflow.Powerflow) ->
         bool

      Check to see if there are any Q limit violations in the network

      1. checkQlimViolations(self: gridpack.powerflow.Powerflow, arg0:
         int) -> bool

      Check to see if there are any Q limit violations in an area

   **checkVoltageViolations(*args, **kwargs)**

      Overloaded function.

      1. checkVoltageViolations(self: gridpack.powerflow.Powerflow) ->
         bool

      Check to see if there are any voltage violations in the network”

      Returns:
         True if no violations found

      1. checkVoltageViolations(self: gridpack.powerflow.Powerflow,
         arg0: int) -> bool

      Check to see if there are any voltage violations in an area

      Parameters:
         area (int): area index

      Returns:
         True if no violations found

   **clearQlimViolations(self: gridpack.powerflow.Powerflow) -> None**

   **clearVoltageViolations(self: gridpack.powerflow.Powerflow) ->
   None**

      Clear “ignore” parameter on all buses

   **close(self: gridpack.powerflow.Powerflow) -> None**

      Close redirect output file

   **exportPSSE23(self: gridpack.powerflow.Powerflow, arg0: str) ->
   None**

      Export final configuration to PSS/E v23 formatted file

   **exportPSSE33(self: gridpack.powerflow.Powerflow, arg0: str) ->
   None**

      Export final configuration to PSS/E v33 formatted file

   **exportPSSE34(self: gridpack.powerflow.Powerflow, arg0: str) ->
   None**

      Export final configuration to PSS/E v34 formatted file

   **getDataCollectionBranchParamInt(self:
   gridpack.powerflow.Powerflow, arg0: int, arg1: int, arg2: str,
   arg3: str) -> object**

      Get (integer) load parameters in data collection for specified
      bus

   **getDataCollectionBranchParamReal(self:
   gridpack.powerflow.Powerflow, arg0: int, arg1: int, arg2: str,
   arg3: str) -> object**

      Get (real) load parameters in data collection for specified bus

   **getDataCollectionBusParamInt(self: gridpack.powerflow.Powerflow,
   arg0: int, arg1: str) -> object**

      Get (integer) load parameters in data collection for specified
      bus

   **getDataCollectionBusParamReal(self: gridpack.powerflow.Powerflow,
   arg0: int, arg1: str) -> object**

      Get (real) load parameters in data collection for specified bus

   **getDataCollectionGenParamInt(self: gridpack.powerflow.Powerflow,
   arg0: int, arg1: str, arg2: str) -> object**

      Get (integer) generator parameters in data collection for
      specified bus

   **getDataCollectionGenParamReal(self: gridpack.powerflow.Powerflow,
   arg0: int, arg1: str, arg2: str) -> object**

      Get (real) generator parameters in data collection for specified
      bus

   **getDataCollectionLoadParamInt(self: gridpack.powerflow.Powerflow,
   arg0: int, arg1: str, arg2: str) -> object**

      Get (integer) load parameters in data collection for specified
      bus

   **getDataCollectionLoadParamReal(self:
   gridpack.powerflow.Powerflow, arg0: int, arg1: str, arg2: str) ->
   object**

      Get (real) load parameters in data collection for specified bus

   **getGeneratorMargins(self: gridpack.powerflow.Powerflow, arg0:
   int, arg1: int) -> object**

      Get the current real power generation and the maximum and
      minimum total power generation for all generators in a zone. If
      zone is less than 1 then return values for all generators in the
      area.

   **getPFSolutionSingleBus(self: gridpack.powerflow.Powerflow, arg0:
   int) -> object**

      get the power flow solution for the specific bus, vmag and v
      angle

   **getTotalLoadRealPower(self: gridpack.powerflow.Powerflow, arg0:
   int, arg1: int) -> float**

   **ignoreVoltageViolations(self: gridpack.powerflow.Powerflow) ->
   None**

      Set “ignore” parameter on all buses with violations so that
      subsequent checks are not counted as violations

   **initialize(self: gridpack.powerflow.Powerflow) -> None**

      Initialize the power network. ``readNetwork()`` must be  called
      before this.

   **modifyDataCollectionBranchParam(*args, **kwargs)**

      Overloaded function.

      1. modifyDataCollectionBranchParam(self:
         gridpack.powerflow.Powerflow, arg0: int, arg1: int, arg2:
         str, arg3: str, arg4: float) -> bool

      Modify (real) parameters in data collection for specified branch

      1. modifyDataCollectionBranchParam(self:
         gridpack.powerflow.Powerflow, arg0: int, arg1: int, arg2:
         str, arg3: str, arg4: int) -> bool

      Modify (integer) parameters in data collection for specified
      branch

   **modifyDataCollectionBusParam(*args, **kwargs)**

      Overloaded function.

      1. modifyDataCollectionBusParam(self:
         gridpack.powerflow.Powerflow, arg0: int, arg1: str, arg2:
         float) -> bool

      Modify (real) parameters in data collection for specified bus

      1. modifyDataCollectionBusParam(self:
         gridpack.powerflow.Powerflow, arg0: int, arg1: str, arg2:
         int) -> bool

      Modify (integer) parameters in data collection for specified bus

   **modifyDataCollectionGenParam(*args, **kwargs)**

      Overloaded function.

      1. modifyDataCollectionGenParam(self:
         gridpack.powerflow.Powerflow, arg0: int, arg1: str, arg2:
         str, arg3: float) -> bool

      Modify (real) generator parameters in data collection for
      specified bus

      1. modifyDataCollectionGenParam(self:
         gridpack.powerflow.Powerflow, arg0: int, arg1: str, arg2:
         str, arg3: int) -> bool

      Modify (integer) generator parameters in data collection for
      specified bus

   **modifyDataCollectionLoadParam(*args, **kwargs)**

      Overloaded function.

      1. modifyDataCollectionLoadParam(self:
         gridpack.powerflow.Powerflow, arg0: int, arg1: str, arg2:
         str, arg3: float) -> bool

      Modify (real) load parameters in data collection for specified
      bus

      1. modifyDataCollectionLoadParam(self:
         gridpack.powerflow.Powerflow, arg0: int, arg1: str, arg2:
         str, arg3: int) -> bool

      Modify (integer) load parameters in data collection for
      specified bus

   **nl_solve(self: gridpack.powerflow.Powerflow) -> bool**

      Solve the power flow problem using the math library non-linear
      solver

   **open(self: gridpack.powerflow.Powerflow, arg0: str) -> None**

      Redirect power flow module output from standard out to a file

   **print(self: gridpack.powerflow.Powerflow, arg0: str) -> None**

      Print an arbitrary string to output

   **readNetwork(self: gridpack.powerflow.Powerflow, arg0:
   gridpack.Configuration, arg1: int) -> None**

      Read the network specified in the configuration.

      Parameters:
         config (gridpack.Configuration): power flow problem
         configuration usually read from an XML file idx (int): The
         index of the power flow problem in the configuration

   **reload(self: gridpack.powerflow.Powerflow) -> None**

      Reload network state and parameters from data collections

   **resetPower(self: gridpack.powerflow.Powerflow) -> None**

      Reset power of loads and generators to original values

   **resetVoltages(self: gridpack.powerflow.Powerflow) -> None**

   **saveData(self: gridpack.powerflow.Powerflow) -> None**

      Save results of powerflow calculation to network data collection
      objects

   **saveDataAlsotoOrg(self: gridpack.powerflow.Powerflow) -> None**

      Save results of powerflow calculation to data collection objects
      also modify the original bus mag, ang, and the original
      generator PG QG in the datacollection

   **scaleGeneratorRealPower(self: gridpack.powerflow.Powerflow, arg0:
   float, arg1: int, arg2: int) -> None**

      Scale load power.

      Parameters:
         scale (float): factor to scale real power generation area
         (int): area index of area for scaling generation zone (int):
         zone index of zone for scaling generation

   **scaleLoadPower(self: gridpack.powerflow.Powerflow, arg0: float,
   arg1: int, arg2: int) -> None**

   **setVoltageLimits(self: gridpack.powerflow.Powerflow, arg0: float,
   arg1: float) -> None**

      “Set voltage limits on all buses”

      Parameters:
         Vmin (float): lower bound on voltages Vmax (float): upper
         bound on voltages

   **solve(self: gridpack.powerflow.Powerflow) -> bool**

      Solve the power flow problem using a custom Newton-Raphson loop

   **suppressOutput(self: gridpack.powerflow.Powerflow, arg0: bool) ->
   None**

      Suppress all output from power flow module

   **useRateB(self: gridpack.powerflow.Powerflow, arg0: bool) ->
   None**

      Use rate B parameter for line overload violations

   **write(self: gridpack.powerflow.Powerflow) -> None**

      Write out results of powerflow calculation to standard output

   **writeBranch(self: gridpack.powerflow.Powerflow, arg0: str) ->
   None**

      Write branch results

   **writeBus(self: gridpack.powerflow.Powerflow, arg0: str) -> None**

      Write bus results

   **writeHeader(self: gridpack.powerflow.Powerflow, arg0: str) ->
   None**

      Write the specified header string to output

   **writeRTPRDiagnostics(self: gridpack.powerflow.Powerflow, arg0:
   int, arg1: int, arg2: int, arg3: int, arg4: int, arg5: int, arg6:
   str) -> None**

      Write real time path rating diagnostics
