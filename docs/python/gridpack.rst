======
Base
======

The base ``gridpack`` Python module contains several utility classes.


Parallel Environment
======================

The GridPACK parallel environment is initialized by instantiating an
instance of ``gridpack.Environment``.  

Usage
-----

.. code:: python
   
   import gridpack
   
   env = gridpack.Environment()
   comm = gridpack.Communicator()

   c = gridpack.Communicator()
   sys.stdout.write("hello from process %d of %d\n" %
                    (c.rank(), c.size()))


Reference
---------

.. autoclass:: gridpack::Environment
   :members:

.. autoclass:: gridpack::Communicator
   :members:
             

Utility Classes
===============

Reference
---------

.. autoclass:: gridpack::Configuration
   :members:
.. autoclass:: gridpack::CoarseTimer
   :members:
.. autoclass:: gridpack::NoPrint
   :members:


Task Manager
============

Usage
-----

Reference
---------

.. autoclass:: gridpack::TaskCounter
   :members:

.. autoclass:: gridpack::TaskManager
   :members:
