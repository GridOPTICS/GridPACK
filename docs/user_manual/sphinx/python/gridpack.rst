
Base
****

The base ``gridpack`` Python module contains several utility classes.


Parallel Environment
====================

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

**class gridpack.Environment**

   GridPACK parallel environment

**class gridpack.Communicator**

   GridPACK parallel communicator

   **barrier(self: gridpack.Communicator) -> None**

   **divide(self: gridpack.Communicator, arg0: int) ->
   `gridpack.Communicator <#gridpack.Communicator>`_**

   **rank(self: gridpack.Communicator) -> int**

   **size(self: gridpack.Communicator) -> int**

   **split(self: gridpack.Communicator, arg0: int) ->
   `gridpack.Communicator <#gridpack.Communicator>`_**

   **sync(self: gridpack.Communicator) -> None**

   **worldRank(self: gridpack.Communicator) -> int**


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

**class gridpack.Configuration**

   Configuration database

   ``KeySep = '.'``

   **get(self: gridpack.Configuration, arg0: str) -> object**

   **getCursor(self: gridpack.Configuration, arg0: str) ->
   gridpack.ConfigurationCursor**

   **open(self: gridpack.Configuration, arg0: str, arg1:
   gridpack.Communicator) -> bool**


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

**class gridpack.TaskCounter**

   ``property task_id``

**class gridpack.TaskManager**

   **cancel(self: gridpack.TaskManager) -> None**

   **nextTask(*args, **kwargs)**

      Overloaded function.

      1. nextTask(self: gridpack.TaskManager, arg0:
         gridpack.TaskCounter) -> bool

      Get the next task from the task manager.

      1. nextTask(self: gridpack.TaskManager, arg0:
         gridpack.Communicator, arg1: gridpack.TaskCounter) -> bool

      Get the next task from the task manager on a (sub)communicator.

   **printStats(self: gridpack.TaskManager) -> None**

      Print out statistics on how tasks are distributed on processors

   **set(self: gridpack.TaskManager, arg0: int) -> None**

      Specify total number of tasks and set task manager to zero


Utility Classes
===============


Reference
---------

**class gridpack.CoarseTimer**

   A general purpose execution timer

   **configTimer(self: gridpack.CoarseTimer, arg0: bool) -> None**

   **createCategory(self: gridpack.CoarseTimer, arg0: str) -> int**

   **currentTime(self: gridpack.CoarseTimer) -> float**

   **dump(self: gridpack.CoarseTimer) -> None**

   **dumpProfile(*args, **kwargs)**

      Overloaded function.

      1. dumpProfile(self: gridpack.CoarseTimer, arg0: int) -> None

      2. dumpProfile(self: gridpack.CoarseTimer, arg0: str) -> None

   **start(self: gridpack.CoarseTimer, arg0: int) -> None**

   **stop(self: gridpack.CoarseTimer, arg0: int) -> None**

**class gridpack.NoPrint**

   Device to control output volume

   **setStatus(self: gridpack.NoPrint, arg0: bool) -> None**

      Set status

   **status(self: gridpack.NoPrint) -> bool**

      Get current status
