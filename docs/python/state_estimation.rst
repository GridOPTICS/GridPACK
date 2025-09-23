==================
 State Estimation
==================

The state estimation Python module closely mirrors the ``SEAppModule``
C++ interface, including reading, manipulation, and setting
``Measurement``s. 

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

.. automodule:: gridpack.state_estimation
   :members:
   :undoc-members:
