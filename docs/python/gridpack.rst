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

A minimal example:

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


Configuration
=============

``Configuration`` is a hierarchical database for storing and
retrieving keyword/value pairs. It is by GridPACK applications to read
initial input files from XML format files.

Usage
-----

A minimal example


Reference
---------

.. autoclass:: gridpack.Configuration
   :members:

Task Manager
============

The ``TaskManager`` is used to distribute and execute an arbitrary
numbre of tasks across available processors.  


Usage
-----

A minimal example

.. code:: python

   import gridpack
   
   env = gridpack.Environment()
   c = gridpack.Communicator()
   tskmgr = gridpack.TaskManager(c)
   
   task = gridpack.TaskCounter()

   tskmgr.set(100)
   while tskmgr.nextTask(task):
       sys.stdout.write("process %d of %d executing task %d\n" %
                        (c.rank(), c.size(), task.task_id))
   tskmgr = None


Reference
---------

.. autoclass:: gridpack::TaskCounter
   :members:

.. autoclass:: gridpack::TaskManager
   :members:


Utility Classes
===============

Reference
---------

.. autoclass:: gridpack::CoarseTimer
   :members:
.. autoclass:: gridpack::NoPrint
   :members:


