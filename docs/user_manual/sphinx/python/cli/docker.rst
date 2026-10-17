Docker
======

Building GridPACK and its dependencies is the hard part of using it.  The
published ``pnnl/gridpack`` image removes that step: it holds GridPACK,
PETSc, Global Arrays, Boost and MPI, and its ``Dockerfile`` installs the
``python/`` package, so any image built from a release that includes this
package also carries the ``gridpack`` command and the Python API::

    docker pull pnnl/gridpack:latest

For most users this is the recommended way to run both the CLI and the
Python API.  Nothing in the package knows it is in a container; the same
commands work on a native install (:doc:`installation`).  This page uses
``pygridpack`` as the image name; substitute ``pnnl/gridpack:latest`` if you
pulled rather than built.

Two ways to get an image with the package in it:

- **Pull** ``pnnl/gridpack:latest`` (or a version tag) from a release made
  after the package was merged.  Check with ``gridpack --version`` inside
  the container; an image from before the merge has an older ``gridpack``
  module (version 0.0.1) and no ``gridpack`` command.
- **Build** ``Dockerfile.pygridpack`` from a checkout, described next.  It
  layers the checkout's ``python/`` onto the published image in minutes, and
  is how to try a branch before it is released.

What the image contains
-----------------------

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - What
     - Where
   * - The ``gridpack`` command and the ``gridpack`` Python package
     - on ``PATH`` and importable by ``python3``, no environment variables needed
   * - Python and MPI
     - Python 3.13, Open MPI 5 (``mpiexec``)
   * - GridPACK C++ drivers (``pf.x``, ``ca.x``, ``dsf.x``, ``dsf2.x``, ``hadrec.x``, ``stes.x``, ``emt.x``)
     - ``/app/src/install/bin``
   * - The example networks and inputs from the source tree
     - ``/app/src/applications/data_sets`` (``input/``, ``raw/``, ``dyr/``, ``contingencies/``, ``emt/``)
   * - The test suite
     - ``/app/python/tests``
   * - Your files
     - ``/app/workspace``, the directory you mount when starting the container

The demos under ``python/examples/`` are not baked in; run them from a
mounted checkout (see below).

Build the image from a checkout
-------------------------------

From the root of a GridPACK checkout::

    git submodule update --init python/pybind11
    docker build -f Dockerfile.pygridpack -t pygridpack .

The first build downloads the base image (several gigabytes); after that the
Python layer builds in minutes.  The build refuses to produce an image whose
import, ``gridpack --version`` or power-flow smoke run fails, then runs the
whole test suite with ``--fail-on-skip``, so a finished image is a tested
one.

To hand a built image to someone without a checkout::

    docker save pygridpack | gzip > pygridpack.tar.gz     # sender
    docker load < pygridpack.tar.gz                        # receiver

Start a container
-----------------

Mount the directory holding your inputs at ``/app/workspace``; it is the
container's working directory, so files you put there are visible to
GridPACK and files GridPACK writes appear back on your machine::

    docker run --rm -it -v "$PWD:/app/workspace" pygridpack

On Windows use PowerShell (``-v "${PWD}:/app/workspace"``) or a WSL2 shell;
inside the container everything is Linux and the rest of this page applies
unchanged.

The container runs as ``root`` by default, so on a Linux host the files it
writes into the mount are root-owned.  Run as yourself instead::

    docker run --rm -it --user $(id -u):$(id -g) -v "$PWD:/app/workspace" pygridpack

To run one command and exit rather than open a shell, append it::

    docker run --rm -v "$PWD:/app/workspace" pygridpack gridpack pf input_14.xml

Prepare a case
--------------

Every GridPACK application reads one XML file, and the XML names its network
(RAW) and dynamics (DYR) files by bare filename, resolved against the working
directory.  Copy the XML and the files it names into one directory and run
there.  A missing network file fails with ``wrong dimension specified``
rather than "file not found".  Starting from the shipped IEEE-14 case::

    D=/app/src/applications/data_sets
    mkdir ieee14 && cd ieee14
    cp $D/input/powerflow/input_14.xml $D/raw/IEEE14.raw .

The EMT application also needs the PETSc options file that ships next to
its inputs, or it aborts inside PETSc::

    cp $D/emt/9-bus/{input_9bus2.xml,case9_PSCAD.raw,case9_GAST.dyr} $D/emt/.petscrc .

Run the command line
--------------------

::

    gridpack pf input_14.xml                      # power flow
    gridpack pf input_14.xml -q --no-timer        # quiet, no timing report
    mpiexec -np 4 gridpack pf input_14.xml        # four ranks

The subcommands are ``powerflow`` (``pf``), ``dsf``, ``se``, ``hadrec``,
``emt`` and ``ca`` (``contingency``); :doc:`commands` documents each.  The
exit code is 0 for a converged run, 1 for a runtime error, 2 for a usage or
configuration error, and 3 when the solver did not converge, so the command
works in scripts::

    gridpack pf input_14.xml -q --no-timer || echo "power flow failed: $?"

``mpiexec`` inside the container sets the rank count; Docker's own replica
or CPU options do not.  Keep a few buses of work per rank: a network too
small for the rank count hangs rather than failing (the 9-bus cases run on
up to three ranks; IEEE 14 and 118 are comfortable at four).

Run Python
----------

Write the script into the mounted directory and run it from there, serially
or under ``mpiexec``::

    python3 my_study.py
    mpiexec -np 2 python3 my_study.py

A first script, in the ``ieee14`` directory prepared above:

.. code:: python

   import gridpack

   with gridpack.Session() as s:                       # one per process
       pf = gridpack.PowerFlow(s, "input_14.xml", suppress_output=True)
       result = pf.solve()
       for bus in result.buses():                      # collective: every rank
           if s.rank == 0:
               print(bus["busId"], bus["name"], bus["voltage"], bus["angle"])
       if s.rank == 0:
           result.to_csv("buses.csv")

``buses.csv`` lands in your mounted directory.  The :doc:`../getting_started`
page walks through every application from here.

Run the demos and the drivers
-----------------------------

With a checkout mounted, the five demos run in place::

    cd /app/workspace/python/examples
    python3 01_power_flow.py
    mpiexec -np 4 python3 02_contingency_screening.py

Each stages its inputs under ``$TMPDIR`` and prints the path of what it
writes; ``python/examples/README.md`` lists the verified rank count per demo.

The C++ drivers are in the image too, which makes it easy to check the
package against them on the same inputs::

    /app/src/install/bin/pf.x input_14.xml
    gridpack pf input_14.xml --no-timer

The result tables should agree; the driver additionally prints a network
properties block and timing.  For dynamic simulation compare against
``dsf2.x``, whose integration loop the package uses; ``dsf.x`` follows a
different loop and differs in the fifth significant digit.

Run the tests
-------------

The suite runs against the inputs shipped in the image::

    docker run --rm pygridpack pytest /app/python/tests -q -rs --fail-on-skip

A skipped test fails this command on purpose: here a skip means ``mpiexec``
is missing, a test dependency is absent or the extension did not build.

Troubleshooting
---------------

``wrong dimension specified``
    The RAW file named in the XML is not in the working directory.

``Error: no contingencies to run``
    The ``Contingency_analysis`` block names no ``<contingencyList>`` and
    sets neither ``FullBranchN1`` nor ``FullGeneratorN1``.

A run that never finishes, with every rank at full CPU
    Too many ranks for the network.  Reduce ``-np``.

``gridpack emt`` or ``emt.x`` aborts inside PETSc
    Copy ``.petscrc`` from ``/app/src/applications/data_sets/emt/`` into the
    working directory.

Files in the mount owned by root
    Start the container with ``--user $(id -u):$(id -g)``.
