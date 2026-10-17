Python API Reference
====================

The classes below are the supported way to drive GridPACK from Python;
:doc:`getting_started` shows them in use.  Each wraps one GridPACK
application and handles MPI, process shutdown and result gathering.  The
pages under *Low-level bindings* document the pybind11 objects underneath,
which stay available for anything these classes do not cover.

Calls marked *collective* must be made on every rank, or the run deadlocks.

Session
-------

.. autoclass:: gridpack.Session
   :members:

Power flow
----------

.. autoclass:: gridpack.PowerFlow
   :members:

.. autoclass:: gridpack.PowerFlowResult
   :members:

Contingency analysis
--------------------

.. autoclass:: gridpack.ContingencyAnalysis
   :members:

.. autoclass:: gridpack.Contingency
   :members:

.. autoclass:: gridpack.ContingencyResult
   :members:

.. autofunction:: gridpack.contingency.parse_contingency_list

.. autofunction:: gridpack.contingency.generate_n1_contingencies

Dynamic simulation
------------------

.. autoclass:: gridpack.DynamicSim
   :members:

.. autoclass:: gridpack.DynamicSimStepper
   :members:

.. autoclass:: gridpack.DSFResult
   :members:

HADREC
------

.. autoclass:: gridpack.Hadrec
   :members:

.. autofunction:: gridpack.hadrec.make_line_trip_action

.. autofunction:: gridpack.hadrec.make_generator_trip_action

.. autofunction:: gridpack.hadrec.make_load_shed_action

State estimation
----------------

.. autoclass:: gridpack.StateEstimation
   :members:

.. autoclass:: gridpack.StateEstimationResult
   :members:

.. autofunction:: gridpack.state_estimation.read_measurements

Errors
------

.. autoexception:: gridpack.GridPACKError

.. autoexception:: gridpack.PowerFlowDiverged
