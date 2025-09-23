
State Estimation
****************

The state estimation Python module closely mirrors the ``SEAppModule``
C++ interface, including reading, manipulation, and setting
``Measurement`` instances.


Usage
=====

A minimal example

.. code:: python

   import sys, os
   from optparse import OptionParser
   import gridpack
   import gridpack.state_estimation

   env = gridpack.Environment()
   comm = gridpack.Communicator()

   config = gridpack.Configuration()
   config.open(inname, comm)

   se = gridpack.state_estimation.SEApp()

   se.readNetwork(config)
   se.initialize()
   se.readMeasurements()
   se.solve()
   se.saveData()
   se.write()

   del se
   del config
   del comm
   del env


Reference
=========

GridPACK state estimation module

**class gridpack.state_estimation.Measurement**

   Container for state estimation measurement data

   ``property busid``

   **ckt(*args, **kwargs)**

      Overloaded function.

      1. ckt(self: gridpack.state_estimation.Measurement) -> str

      Get measurement ckt.

      1. ckt(self: gridpack.state_estimation.Measurement, arg0: str)
         -> None

      Set measurment ckt”

   ``property deviation``

   ``property fbusid``

   ``property tbusid``

   **type(*args, **kwargs)**

      Overloaded function.

      1. type(self: gridpack.state_estimation.Measurement) -> str

      Get measurement type.

      1. type(self: gridpack.state_estimation.Measurement, arg0: str)
         -> None

      Set measurment type, one of “VA”, “VM”, “VUL”, “VLL”, “PIJ”,
      “QIJ”

   ``property value``

**class gridpack.state_estimation.MeasurementVector**

   **append(self: gridpack.state_estimation.MeasurementVector, x:
   gridpack.state_estimation.Measurement) -> None**

      Add an item to the end of the list

   **clear(self: gridpack.state_estimation.MeasurementVector) ->
   None**

      Clear the contents

   **extend(*args, **kwargs)**

      Overloaded function.

      1. extend(self: gridpack.state_estimation.MeasurementVector, L:
         gridpack.state_estimation.MeasurementVector) -> None

      Extend the list by appending all the items in the given list

      1. extend(self: gridpack.state_estimation.MeasurementVector, L:
         Iterable) -> None

      Extend the list by appending all the items in the given list

   **insert(self: gridpack.state_estimation.MeasurementVector, i: int,
   x: gridpack.state_estimation.Measurement) -> None**

      Insert an item at a given position.

   **pop(*args, **kwargs)**

      Overloaded function.

      1. pop(self: gridpack.state_estimation.MeasurementVector) ->
         gridpack.state_estimation.Measurement

      Remove and return the last item

      1. pop(self: gridpack.state_estimation.MeasurementVector, i:
         int) -> gridpack.state_estimation.Measurement

      Remove and return the item at index ``i``

**class gridpack.state_estimation.SEApp**

   State estimation application module

   **addVoltageLimitMeasurements(self:
   gridpack.state_estimation.SEApp, vmin: float = 0.9, vmax: float =
   1.1, deviation: float = 0.001) -> None**

      Add virtual measurements to enforce voltage magnitude
      constraints

      Parameters:
         vmin (float): minimum allowed voltage magnitude vmax (float):
         maximum allowed voltage magnitude deviation (float): standard
         deviation for virtual measurements

   **checkMeasurementConsistency(self:
   gridpack.state_estimation.SEApp) -> None**

      Check for potential measurement inconsistencies Identifies cases
      where measurements may conflict with physical constraints

   **debugMapper(self: gridpack.state_estimation.SEApp) -> None**

      Perform a targeted check for bad measurements on Bus 8

   **getMeasurements(self: gridpack.state_estimation.SEApp, arg0:
   gridpack.Configuration) ->
   `gridpack.state_estimation.MeasurementVector
   <#gridpack.state_estimation.MeasurementVector>`_**

      Get measurements from an open Configuration.

   **handlePVBusVoltages(self: gridpack.state_estimation.SEApp) ->
   None**

      Apply proper treatment for PV bus voltage measurements Ensures
      PV bus voltages are treated as constraints rather than regular
      measurements

   **handleVAMeasurements(self: gridpack.state_estimation.SEApp) ->
   None**

      Apply special handling for voltage angle (VA) measurements
      Ensures angle measurements at slack buses are treated as
      constraints

   **hasConverged(self: gridpack.state_estimation.SEApp) -> bool**

   **identifyPVBusConstraints(self: gridpack.state_estimation.SEApp)
   -> None**

      Identify PV buses in the network and their connections Used to
      properly handle voltage constraints at generator buses

   **initialize(self: gridpack.state_estimation.SEApp) -> None**

      Set up exchange buffers and other internal parameters and
      initialize network components using data from data collection

   **preCheckMeasurements(self: gridpack.state_estimation.SEApp) ->
   None**

      Perform pre-check for suspicious measurements

   **readMeasurements(self: gridpack.state_estimation.SEApp) -> None**

      Read branch and bus measurements. These will come from a
      separate file.  The name of this file comes from the input
      configuration from which the network was read.

   **readNetwork(self: gridpack.state_estimation.SEApp, arg0:
   gridpack.Configuration) -> None**

      Read in and partition the network. The input file is read
      directly from the state_estimation block in the specified
      configuration.

      Parameters:
         config (gridpack.Configuration): power flow problem
         configuration usually read from an XML file

   **reportJacobianPerformance(self: gridpack.state_estimation.SEApp)
   -> None**

      Report Jacobian optimization performance statistics

   **saveData(self: gridpack.state_estimation.SEApp) -> None**

      Save results of state estimation calculation to data collection
      objects

   **setMeasurements(self: gridpack.state_estimation.SEApp, arg0:
   gridpack.state_estimation.MeasurementVector) -> None**

      Set bus and branch measurements.

      Parameters:

      measurments a list of Measurement possibly from
      getMeasurements()

   **solve(self: gridpack.state_estimation.SEApp) -> None**

      Solve the state estimation problem

   **write(self: gridpack.state_estimation.SEApp) -> None**

      Write final results of state estimation calculation to standard
      output
