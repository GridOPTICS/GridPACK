
HADREC
******


Reference
=========

GridPACK HADREC Application module

**class gridpack.hadrec.Action**

   ``property actiontype``

   ``property branch_ckt``

   ``property brch_from_bus_number``

   ``property brch_to_bus_number``

   ``property bus_number``

   ``property componentID``

   ``property percentage``

**class gridpack.hadrec.Module**

   **applyAction(self: gridpack.hadrec.Module, arg0:
   gridpack.hadrec.Action) -> None**

   **executeDynSimuOneStep(self: gridpack.hadrec.Module) -> None**

   **exportPSSE23(self: gridpack.hadrec.Module, s: str = '') -> None**

   **exportPSSE33(self: gridpack.hadrec.Module, s: str = '') -> None**

   **exportPSSE34(self: gridpack.hadrec.Module, s: str = '') -> None**

   **fullInitializationBeforeDynSimuSteps(self:
   gridpack.hadrec.Module, s: str = '', BusFaults:
   gridpack.dynamic_simulation.EventVector =
   <gridpack.dynamic_simulation.EventVector object at 0x1296345e2eb0>,
   pfcase_idx: int = -1, dscase_idx: int) -> None**

   **getBranchEndpoints(self: gridpack.hadrec.Module, arg0: int) ->
   tuple**

   **getBranchInfoBool(self: gridpack.hadrec.Module, bus_idx: int,
   name: str, dev_idx: int | None = -1) -> object**

   **getBranchInfoInt(self: gridpack.hadrec.Module, bus_idx: int,
   name: str, dev_idx: int | None = -1) -> object**

   **getBranchInfoReal(self: gridpack.hadrec.Module, bus_idx: int,
   name: str, dev_idx: int | None = -1) -> object**

   **getBranchInfoString(self: gridpack.hadrec.Module, bus_idx: int,
   name: str, dev_idx: int | None = -1) -> object**

   **getBusInfoBool(self: gridpack.hadrec.Module, bus_idx: int, name:
   str, dev_idx: int | None = -1) -> object**

   **getBusInfoInt(self: gridpack.hadrec.Module, bus_idx: int, name:
   str, dev_idx: int | None = -1) -> object**

   **getBusInfoReal(self: gridpack.hadrec.Module, bus_idx: int, name:
   str, dev_idx: int | None = -1) -> object**

   **getBusInfoString(self: gridpack.hadrec.Module, bus_idx: int,
   name: str, dev_idx: int | None = -1) -> object**

   **getBusTotalLoadPower(self: gridpack.hadrec.Module, arg0: int) ->
   object**

   **getConnectedBranches(self: gridpack.hadrec.Module, arg0: int) ->
   list[int]**

   **getDataCollectionBranchParam(*args, **kwargs)**

      Overloaded function.

      1. getDataCollectionBranchParam(self: gridpack.hadrec.Module,
         arg0: int, arg1: int, arg2: str, arg3: str) -> object

      2. getDataCollectionBranchParam(self: gridpack.hadrec.Module,
         arg0: int, arg1: int, arg2: str, arg3: str) -> object

   **getDataCollectionBusParam(*args, **kwargs)**

      Overloaded function.

      1. getDataCollectionBusParam(self: gridpack.hadrec.Module, arg0:
         int, arg1: str) -> object

      2. getDataCollectionBusParam(self: gridpack.hadrec.Module, arg0:
         int, arg1: str) -> object

   **getDataCollectionGenParam(*args, **kwargs)**

      Overloaded function.

      1. getDataCollectionGenParam(self: gridpack.hadrec.Module, arg0:
         int, arg1: str, arg2: str) -> object

      2. getDataCollectionGenParam(self: gridpack.hadrec.Module, arg0:
         int, arg1: str, arg2: str) -> object

   **getDataCollectionLoadParam(*args, **kwargs)**

      Overloaded function.

      1. getDataCollectionLoadParam(self: gridpack.hadrec.Module,
         arg0: int, arg1: str, arg2: str) -> object

      2. getDataCollectionLoadParam(self: gridpack.hadrec.Module,
         arg0: int, arg1: str, arg2: str) -> object

   **getGeneratorPower(self: gridpack.hadrec.Module, arg0: int, arg1:
   str) -> object**

   **getObservationLists(self: gridpack.hadrec.Module) -> tuple**

   **getObservationLists_withBusFreq(self: gridpack.hadrec.Module) ->
   tuple**

   **getObservations(self: gridpack.hadrec.Module) -> list[float]**

   **getPFSolutionSingleBus(self: gridpack.hadrec.Module, arg0: int)
   -> object**

   **getState(self: gridpack.hadrec.Module, arg0: int, arg1: str,
   arg2: str, arg3: str) -> object**

   **getZoneGeneratorPower(self: gridpack.hadrec.Module) -> object**

   **getZoneLoads(self: gridpack.hadrec.Module) -> object**

   **initializeDynSimu(self: gridpack.hadrec.Module, faults:
   gridpack.dynamic_simulation.EventVector =
   <gridpack.dynamic_simulation.EventVector object at 0x1296345cc870>,
   dscase_idx: int = -1) -> None**

   **isDynSimuDone(self: gridpack.hadrec.Module) -> bool**

   **modifyDataCollectionBranchParam(*args, **kwargs)**

      Overloaded function.

      1. modifyDataCollectionBranchParam(self: gridpack.hadrec.Module,
         arg0: int, arg1: int, arg2: str, arg3: str, arg4: float) ->
         bool

      2. modifyDataCollectionBranchParam(self: gridpack.hadrec.Module,
         arg0: int, arg1: int, arg2: str, arg3: str, arg4: int) ->
         bool

   **modifyDataCollectionBusParam(*args, **kwargs)**

      Overloaded function.

      1. modifyDataCollectionBusParam(self: gridpack.hadrec.Module,
         arg0: int, arg1: str, arg2: float) -> bool

      2. modifyDataCollectionBusParam(self: gridpack.hadrec.Module,
         arg0: int, arg1: str, arg2: int) -> bool

   **modifyDataCollectionGenParam(*args, **kwargs)**

      Overloaded function.

      1. modifyDataCollectionGenParam(self: gridpack.hadrec.Module,
         arg0: int, arg1: str, arg2: str, arg3: float) -> bool

      2. modifyDataCollectionGenParam(self: gridpack.hadrec.Module,
         arg0: int, arg1: str, arg2: str, arg3: int) -> bool

   **modifyDataCollectionLoadParam(*args, **kwargs)**

      Overloaded function.

      1. modifyDataCollectionLoadParam(self: gridpack.hadrec.Module,
         arg0: int, arg1: str, arg2: str, arg3: float) -> bool

      2. modifyDataCollectionLoadParam(self: gridpack.hadrec.Module,
         arg0: int, arg1: str, arg2: str, arg3: int) -> bool

   **numGenerators(*args, **kwargs)**

      Overloaded function.

      1. numGenerators(self: gridpack.hadrec.Module) -> int

      2. numGenerators(self: gridpack.hadrec.Module, arg0: int) -> int

   **numLines(*args, **kwargs)**

      Overloaded function.

      1. numLines(self: gridpack.hadrec.Module) -> int

      2. numLines(self: gridpack.hadrec.Module, arg0: int) -> int

   **numLoads(*args, **kwargs)**

      Overloaded function.

      1. numLoads(self: gridpack.hadrec.Module) -> int

      2. numLoads(self: gridpack.hadrec.Module, arg0: int) -> int

   **numStorage(*args, **kwargs)**

      Overloaded function.

      1. numStorage(self: gridpack.hadrec.Module) -> int

      2. numStorage(self: gridpack.hadrec.Module, arg0: int) -> int

   **readPowerFlowData(self: gridpack.hadrec.Module, s: str = '',
   pfcase_idx: int = -1) -> None**

   **scatterInjectionLoad(self: gridpack.hadrec.Module, vbusNum:
   list[int] = [], vloadP: list[float] = [], vloadQ: list[float] = [])
   -> None**

   **scatterInjectionLoadNew(self: gridpack.hadrec.Module, vbusNum:
   list[int] = [], vloadP: list[float] = [], vloadQ: list[float] = [])
   -> None**

   **scatterInjectionLoadNewConstCur(self: gridpack.hadrec.Module,
   vbusNum: list[int] = [], vCurR: list[float] = [], vCurI:
   list[float] = []) -> None**

   **scatterInjectionLoadNew_Norton(self: gridpack.hadrec.Module,
   vbusNum: list[int] = [], vloadP: list[float] = [], vloadQ:
   list[float] = [], vimpedanceR: list[float] = [], vimpedanceI:
   list[float] = []) -> None**

   **scatterInjectionLoadNew_compensateY(self: gridpack.hadrec.Module,
   vbusNum: list[int] = [], vloadP: list[float] = [], vloadQ:
   list[float] = []) -> None**

   **setState(self: gridpack.hadrec.Module, arg0: int, arg1: str,
   arg2: str, arg3: str, arg4: float) -> bool**

   **setWideAreaControlSignal(self: gridpack.hadrec.Module,
   bus_number: int = '-1', genid: str = '', wideAreaControlSignal:
   float = 0.0) -> None**

   **solvePowerFlow(self: gridpack.hadrec.Module) -> bool**

   **solvePowerFlowBeforeDynSimu(self: gridpack.hadrec.Module, s: str
   = '', pfcase_idx: int = -1) -> None**

   **solvePowerFlowBeforeDynSimu_withFlag(self:
   gridpack.hadrec.Module, s: str = '', pfcase_idx: int = -1) ->
   bool**

   **totalBranches(self: gridpack.hadrec.Module) -> int**

   **totalBuses(self: gridpack.hadrec.Module) -> int**

   **transferPFtoDS(self: gridpack.hadrec.Module) -> None**

   **updateData(self: gridpack.hadrec.Module) -> None**
