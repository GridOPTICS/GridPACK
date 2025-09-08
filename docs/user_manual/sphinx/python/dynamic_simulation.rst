
Dynamic Simulation
******************


Reference
=========

GridPACK Dynamic Simulation Application module

**class gridpack.dynamic_simulation.DSFullApp**

   **applyConstYLoadShedding(self:
   gridpack.dynamic_simulation.DSFullApp, arg0: int, arg1: float) ->
   None**

   **applyConstYLoad_Change_P(self:
   gridpack.dynamic_simulation.DSFullApp, arg0: int, arg1: float) ->
   None**

   **applyConstYLoad_Change_Q(self:
   gridpack.dynamic_simulation.DSFullApp, arg0: int, arg1: float) ->
   None**

   **applyGFIAdjustment(self: gridpack.dynamic_simulation.DSFullApp,
   arg0: int, arg1: int, arg2: str, arg3: float) -> None**

   **applyGeneratorTripping(self:
   gridpack.dynamic_simulation.DSFullApp, arg0: int, arg1: str) ->
   None**

   **applyLoadShedding(self: gridpack.dynamic_simulation.DSFullApp,
   arg0: int, arg1: str, arg2: float) -> None**

   **clearConstYLoad_Change_P(self:
   gridpack.dynamic_simulation.DSFullApp) -> None**

   **clearConstYLoad_Change_Q(self:
   gridpack.dynamic_simulation.DSFullApp) -> None**

   **clearLineTripAction(self: gridpack.dynamic_simulation.DSFullApp)
   -> None**

   **close(self: gridpack.dynamic_simulation.DSFullApp) -> None**

   **executeOneSimuStep(self: gridpack.dynamic_simulation.DSFullApp)
   -> None**

   **frequencyOK(self: gridpack.dynamic_simulation.DSFullApp) ->
   bool**

   **getBranchEndpoints(self: gridpack.dynamic_simulation.DSFullApp,
   arg0: int) -> tuple**

   **getBranchInfoBool(self: gridpack.dynamic_simulation.DSFullApp,
   bus_idx: int, name: str, dev_idx: int | None = -1) -> object**

   **getBranchInfoInt(self: gridpack.dynamic_simulation.DSFullApp,
   bus_idx: int, name: str, dev_idx: int | None = -1) -> object**

   **getBranchInfoReal(self: gridpack.dynamic_simulation.DSFullApp,
   bus_idx: int, name: str, dev_idx: int | None = -1) -> object**

   **getBranchInfoString(self: gridpack.dynamic_simulation.DSFullApp,
   bus_idx: int, name: str, dev_idx: int | None = -1) -> object**

   **getBusInfoBool(self: gridpack.dynamic_simulation.DSFullApp, arg0:
   int, arg1: str, arg2: int | None) -> object**

   **getBusInfoInt(self: gridpack.dynamic_simulation.DSFullApp,
   bus_idx: int, name: str, dev_idx: int | None = -1) -> object**

   **getBusInfoReal(self: gridpack.dynamic_simulation.DSFullApp,
   bus_idx: int, name: str, dev_idx: int | None = -1) -> object**

   **getBusInfoString(self: gridpack.dynamic_simulation.DSFullApp,
   bus_idx: int, name: str, dev_idx: int | None = -1) -> object**

   **getBusTotalLoadPower(self: gridpack.dynamic_simulation.DSFullApp,
   arg0: int) -> object**

   **getConnectedBranches(self: gridpack.dynamic_simulation.DSFullApp,
   arg0: int) -> list[int]**

   **getCurrentTime(self: gridpack.dynamic_simulation.DSFullApp) ->
   float**

   **getEvents(*args, **kwargs)**

      Overloaded function.

      1. getEvents(self: gridpack.dynamic_simulation.DSFullApp) ->
         gridpack.dynamic_simulation.EventVector

      2. getEvents(self: gridpack.dynamic_simulation.DSFullApp, arg0:
         gridpack.ConfigurationCursor) ->
         gridpack.dynamic_simulation.EventVector

   **getFinalTime(self: gridpack.dynamic_simulation.DSFullApp) ->
   float**

   **getFrequencyFailures(self: gridpack.dynamic_simulation.DSFullApp)
   -> list[int]**

   **getGeneratorMargins(self: gridpack.dynamic_simulation.DSFullApp,
   arg0: int, arg1: int) -> object**

   **getGeneratorPower(self: gridpack.dynamic_simulation.DSFullApp,
   arg0: int, arg1: str) -> object**

   **getGeneratorTimeSeries(self:
   gridpack.dynamic_simulation.DSFullApp) -> list[list[float]]**

   **getListWatchedGenerators(self:
   gridpack.dynamic_simulation.DSFullApp) -> object**

   **getObservationLists(self: gridpack.dynamic_simulation.DSFullApp)
   -> object**

   **getObservationLists_withBusFreq(self:
   gridpack.dynamic_simulation.DSFullApp) -> object**

   **getObservations(self: gridpack.dynamic_simulation.DSFullApp) ->
   object**

   **getObservations_withBusFreq(self:
   gridpack.dynamic_simulation.DSFullApp) -> object**

   **getState(self: gridpack.dynamic_simulation.DSFullApp, arg0: int,
   arg1: str, arg2: str, arg3: str) -> object**

   **getTimeSeriesMap(self: gridpack.dynamic_simulation.DSFullApp) ->
   list[int]**

   **getTimeStep(self: gridpack.dynamic_simulation.DSFullApp) ->
   float**

   **getTotalLoadRealPower(self:
   gridpack.dynamic_simulation.DSFullApp, arg0: int, arg1: int) ->
   float**

   **getZoneGeneratorPower(self:
   gridpack.dynamic_simulation.DSFullApp) -> object**

   **getZoneLoads(self: gridpack.dynamic_simulation.DSFullApp) ->
   object**

   **initialize(self: gridpack.dynamic_simulation.DSFullApp) -> None**

   **isDynSimuDone(self: gridpack.dynamic_simulation.DSFullApp) ->
   bool**

   **isSecure(self: gridpack.dynamic_simulation.DSFullApp) -> int**

   **modifyDataCollectionBusParam(*args, **kwargs)**

      Overloaded function.

      1. modifyDataCollectionBusParam(self:
         gridpack.dynamic_simulation.DSFullApp, arg0: int, arg1: str,
         arg2: float) -> bool

      2. modifyDataCollectionBusParam(self:
         gridpack.dynamic_simulation.DSFullApp, arg0: int, arg1: str,
         arg2: int) -> bool

   **modifyDataCollectionGenParam(*args, **kwargs)**

      Overloaded function.

      1. modifyDataCollectionGenParam(self:
         gridpack.dynamic_simulation.DSFullApp, arg0: int, arg1: str,
         arg2: str, arg3: float) -> bool

      2. modifyDataCollectionGenParam(self:
         gridpack.dynamic_simulation.DSFullApp, arg0: int, arg1: str,
         arg2: str, arg3: int) -> bool

   **modifyDataCollectionLoadParam(*args, **kwargs)**

      Overloaded function.

      1. modifyDataCollectionLoadParam(self:
         gridpack.dynamic_simulation.DSFullApp, arg0: int, arg1: str,
         arg2: str, arg3: float) -> bool

      2. modifyDataCollectionLoadParam(self:
         gridpack.dynamic_simulation.DSFullApp, arg0: int, arg1: str,
         arg2: str, arg3: int) -> bool

   **numGenerators(*args, **kwargs)**

      Overloaded function.

      1. numGenerators(self: gridpack.dynamic_simulation.DSFullApp) ->
         int

      2. numGenerators(self: gridpack.dynamic_simulation.DSFullApp,
         arg0: int) -> int

   **numLines(*args, **kwargs)**

      Overloaded function.

      1. numLines(self: gridpack.dynamic_simulation.DSFullApp) -> int

      2. numLines(self: gridpack.dynamic_simulation.DSFullApp, arg0:
         int) -> int

   **numLoads(*args, **kwargs)**

      Overloaded function.

      1. numLoads(self: gridpack.dynamic_simulation.DSFullApp) -> int

      2. numLoads(self: gridpack.dynamic_simulation.DSFullApp, arg0:
         int) -> int

   **numStorage(*args, **kwargs)**

      Overloaded function.

      1. numStorage(self: gridpack.dynamic_simulation.DSFullApp) ->
         int

      2. numStorage(self: gridpack.dynamic_simulation.DSFullApp, arg0:
         int) -> int

   **open(self: gridpack.dynamic_simulation.DSFullApp, arg0: str) ->
   None**

   **print(self: gridpack.dynamic_simulation.DSFullApp, arg0: str) ->
   None**

   **readGenerators(self: gridpack.dynamic_simulation.DSFullApp,
   ds_idx: int = -1) -> None**

   **readSequenceData(self: gridpack.dynamic_simulation.DSFullApp) ->
   None**

   **reload(self: gridpack.dynamic_simulation.DSFullApp) -> None**

   **reset(self: gridpack.dynamic_simulation.DSFullApp) -> None**

   **resetPower(self: gridpack.dynamic_simulation.DSFullApp) -> None**

   **run(*args, **kwargs)**

      Overloaded function.

      1. run(self: gridpack.dynamic_simulation.DSFullApp) -> None

      2. run(self: gridpack.dynamic_simulation.DSFullApp, arg0: float)
         -> None

   **saveTimeSeries(self: gridpack.dynamic_simulation.DSFullApp, arg0:
   bool) -> None**

   **scaleGeneratorRealPower(self:
   gridpack.dynamic_simulation.DSFullApp, arg0: float, arg1: int,
   arg2: int) -> None**

   **scaleLoadPower(self: gridpack.dynamic_simulation.DSFullApp, arg0:
   float, arg1: int, arg2: int) -> None**

   **scatterInjectionLoad(self: gridpack.dynamic_simulation.DSFullApp,
   arg0: list[int], arg1: list[float], arg2: list[float]) -> None**

   **scatterInjectionLoadNew(self:
   gridpack.dynamic_simulation.DSFullApp, arg0: list[int], arg1:
   list[float], arg2: list[float]) -> None**

   **scatterInjectionLoadNewConstCur(self:
   gridpack.dynamic_simulation.DSFullApp, arg0: list[int], arg1:
   list[float], arg2: list[float]) -> None**

   **scatterInjectionLoadNew_Norton(self:
   gridpack.dynamic_simulation.DSFullApp, arg0: list[int], arg1:
   list[float], arg2: list[float], arg3: list[float], arg4:
   list[float]) -> None**

   **scatterInjectionLoadNew_compensateY(self:
   gridpack.dynamic_simulation.DSFullApp, arg0: list[int], arg1:
   list[float], arg2: list[float]) -> None**

   **setConstYLoadImpedance(self:
   gridpack.dynamic_simulation.DSFullApp, arg0: int, arg1: float,
   arg2: float) -> None**

   **setConstYLoadtoZero_P(self:
   gridpack.dynamic_simulation.DSFullApp, arg0: int) -> None**

   **setConstYLoadtoZero_Q(self:
   gridpack.dynamic_simulation.DSFullApp, arg0: int) -> None**

   **setEvent(self: gridpack.dynamic_simulation.DSFullApp, arg0:
   gridpack.dynamic_simulation.Event) -> None**

   **setFinalTime(self: gridpack.dynamic_simulation.DSFullApp, arg0:
   float) -> None**

   **setFrequencyMonitoring(self:
   gridpack.dynamic_simulation.DSFullApp, arg0: bool, arg1: float) ->
   None**

   **setGeneratorWatch(*args, **kwargs)**

      Overloaded function.

      1. setGeneratorWatch(self:
         gridpack.dynamic_simulation.DSFullApp) -> None

      2. setGeneratorWatch(self:
         gridpack.dynamic_simulation.DSFullApp, arg0: str) -> None

      3. setGeneratorWatch(self:
         gridpack.dynamic_simulation.DSFullApp, buses: list[int] = [],
         tags: list[str] = [], writeFile: bool = True) -> None

      4. setGeneratorWatch(self:
         gridpack.dynamic_simulation.DSFullApp, arg0: str, arg1:
         gridpack.ConfigurationCursor) -> None

      5. setGeneratorWatch(self:
         gridpack.dynamic_simulation.DSFullApp, arg0:
         gridpack.ConfigurationCursor) -> None

   **setLineTripAction(*args, **kwargs)**

      Overloaded function.

      1. setLineTripAction(self:
         gridpack.dynamic_simulation.DSFullApp, arg0: int, arg1: int,
         arg2: str) -> None

      2. setLineTripAction(self:
         gridpack.dynamic_simulation.DSFullApp, arg0: int) -> None

   **setLoadWatch(self: gridpack.dynamic_simulation.DSFullApp) ->
   None**

   **setObservations(self: gridpack.dynamic_simulation.DSFullApp,
   arg0: gridpack.ConfigurationCursor) -> None**

   **setState(self: gridpack.dynamic_simulation.DSFullApp, arg0: int,
   arg1: str, arg2: str, arg3: str, arg4: float) -> bool**

   **setTimeStep(self: gridpack.dynamic_simulation.DSFullApp, arg0:
   float) -> None**

   **setWideAreaControlSignal(self:
   gridpack.dynamic_simulation.DSFullApp, arg0: int, arg1: str, arg2:
   float) -> None**

   **setup(self: gridpack.dynamic_simulation.DSFullApp) -> None**

   **solve(self: gridpack.dynamic_simulation.DSFullApp, arg0:
   gridpack.dynamic_simulation.Event) -> None**

   **solvePowerFlowBeforeDynSimu(self:
   gridpack.dynamic_simulation.DSFullApp, arg0: str, arg1: int) ->
   None**

   **solvePreInitialize(self: gridpack.dynamic_simulation.DSFullApp,
   arg0: gridpack.dynamic_simulation.Event) -> None**

   **totalBranches(self: gridpack.dynamic_simulation.DSFullApp) ->
   int**

   **totalBuses(self: gridpack.dynamic_simulation.DSFullApp) -> int**

   **updateData(self: gridpack.dynamic_simulation.DSFullApp) -> None**

   **write(self: gridpack.dynamic_simulation.DSFullApp, arg0: str) ->
   None**

   **writeRTPRDiagnostics(self: gridpack.dynamic_simulation.DSFullApp,
   arg0: int, arg1: int, arg2: int, arg3: int, arg4: float, arg5:
   float, arg6: str) -> None**

**class gridpack.dynamic_simulation.Event**

   ``property bus_idx``

   ``property end``

   ``property from_idx``

   ``property isBus``

   ``property isGenerator``

   ``property isLine``

   ``property start``

   ``property step``

   ``property tag``

   ``property to_idx``

**class gridpack.dynamic_simulation.EventVector**

   **append(self: gridpack.dynamic_simulation.EventVector, x:
   gridpack::dynamic_simulation::Event) -> None**

      Add an item to the end of the list

   **clear(self: gridpack.dynamic_simulation.EventVector) -> None**

      Clear the contents

   **extend(*args, **kwargs)**

      Overloaded function.

      1. extend(self: gridpack.dynamic_simulation.EventVector, L:
         gridpack.dynamic_simulation.EventVector) -> None

      Extend the list by appending all the items in the given list

      1. extend(self: gridpack.dynamic_simulation.EventVector, L:
         Iterable) -> None

      Extend the list by appending all the items in the given list

   **insert(self: gridpack.dynamic_simulation.EventVector, i: int, x:
   gridpack::dynamic_simulation::Event) -> None**

      Insert an item at a given position.

   **pop(*args, **kwargs)**

      Overloaded function.

      1. pop(self: gridpack.dynamic_simulation.EventVector) ->
         gridpack::dynamic_simulation::Event

      Remove and return the last item

      1. pop(self: gridpack.dynamic_simulation.EventVector, i: int) ->
         gridpack::dynamic_simulation::Event

      Remove and return the item at index ``i``
