============
 Power Flow
============

The Python Power Flow interface mirrors closely the  
``PFAppModule`` C++ interface.  See the `that module's documentation
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

   pfapp = gridpack.powerflow.powerflow()
   pfapp.readNetwork(config, -1)
   pfapp.initialize()
   pfapp.solve();
   pfapp.write();

   # try to force deletion order to avoid problems
   del pfapp
   del env


Reference
=========

.. automodule:: gridpack.powerflow
   :members:
   :undoc-members:               
