Getting Started
***************

This page is for someone who has never driven GridPACK from Python.  It
uses the high-level classes in the ``gridpack`` package; the pages that
follow document the pybind11 objects underneath them, which you only need
for things the high-level classes do not expose.


Install
=======

The package is a compiled extension linked against an installed GridPACK,
so GridPACK, PETSc, Global Arrays, Boost and MPI must already be built as
shared libraries.  Then::

    export GRIDPACK_DIR=/path/to/gridpack/install   # holds lib/GridPACK.cmake
    python3 -m venv ~/.venvs/pygridpack
    source ~/.venvs/pygridpack/bin/activate
    pip install setuptools wheel mpi4py
    cd /path/to/GridPACK
    pip install --no-build-isolation -e python/

Check it::

    gridpack --version
    cd python && pytest -q      # every test passes, none skipped

``python/README.md`` in the source tree covers Docker and the details.


Inputs
======

Every application reads one XML file.  The XML names its network (RAW) and
dynamics (DYR) files by bare filename, so copy the XML and the files it
names into one directory and run from there.  The source tree ships cases
under ``src/applications/data_sets/``: XMLs under ``input/``, networks under
``raw/``, machine data under ``dyr/``.  A missing network file fails with
``wrong dimension specified``, not with a file-not-found message.


The Session
===========

Everything starts with one :class:`gridpack.Session` per process.  It
initializes MPI, PETSc and Global Arrays, hands out ``rank`` and ``size``,
and tears down every wrapper you build in the right order when the
``with`` block ends.  Building two in one process raises.

.. code:: python

   import gridpack

   with gridpack.Session() as s:
       print("rank", s.rank, "of", s.size)
       ...   # wrappers go here


Power flow
==========

.. code:: python

   with gridpack.Session() as s:
       pf = gridpack.PowerFlow(s, "input_14.xml", suppress_output=True)
       result = pf.solve()          # raises PowerFlowDiverged if it fails
       print(result.converged, result.iterations)

       for bus in result.buses():   # one dict per bus, every rank
           print(bus["busId"], bus["name"], bus["voltage"], bus["angle"])

       bad = result.violations(min_voltage=0.95, max_voltage=1.05,
                               overload_threshold=100.0)
       print(len(bad["voltage"]), "buses out of band,",
             len(bad["overload"]), "branches overloaded")

       result.to_csv("buses.csv")             # also table="branches" / "generators"

``solve(strict=False)`` returns the unconverged result instead of raising;
``result.mismatch`` then names the worst bus.  ``buses()``, ``branches()``,
``generators()`` and ``violations()`` gather the whole network across ranks
and are **collective**: every rank must call them, so do not put them under
``if s.rank == 0``.


Contingency analysis
====================

Each rank takes contingencies from a shared queue, so more ranks finish the
screen sooner.  The XML must name a ``<contingencyList>`` file; the
``FullBranchN1`` shortcut is not expanded yet.

.. code:: python

   with gridpack.Session() as s:
       with gridpack.ContingencyAnalysis(s, "input_14.xml",
                                         print_calc_files=False,
                                         suppress_output=True) as ca:
           base = ca.solve_base_case()
           ca.run()                        # this rank's share
           results = ca.gather()           # all ranks', sorted by name
       if s.rank == 0:
           for r in results:
               print(r.name.strip(), r.status)

``status`` is one of ``OK``, ``BUS VIOLATION``, ``BRANCH VIOLATION``,
``BUS+BRANCH VIOLATION``, ``DIVERGENT``, ``ISLANDED``, ``SLACK OVERLOAD``
or ``NOT FOUND``, the same verdicts ``ca.x`` prints.


Dynamic simulation
==================

Two routes.  :class:`gridpack.DynamicSim` hands the whole run to C++, with
the fault declared in the XML, and is the path to use when you do not need
to intervene:

.. code:: python

   with gridpack.Session() as s:
       with gridpack.DynamicSim(s, "input_9b3g.xml", suppress_output=True) as ds:
           result = ds.run()
           if s.rank == 0:
               print(result.get_generator_power(1, "1"))   # (P, Q) at the end

The per-step series arrives in the generator watch file the XML names
under ``<generatorWatchFileName>``.

:class:`gridpack.DynamicSimStepper` keeps the loop in Python, so an event
can be decided at run time:

.. code:: python

   with gridpack.Session() as s:
       with gridpack.DynamicSimStepper(s, "input_kundur_two_area_wsieg1.xml",
                                       suppress_output=True) as st:
           while not st.done and st.current_time < 8.0:
               st.step()
               if st.step_count == 500:
                   st.apply_generator_trip(bus_number=2, gen_id="1")
           if s.rank == 0:
               st.result.to_csv("stepped.csv")

The stepper records what the XML's ``<observations>`` block lists;
without one the series is empty.  The XML must also declare at least one
event under ``<Events>``, even one timed past ``simulationTime``.
Channel names look like ``gen_2_1_rspeed`` and ``bus_6_vmag``.


State estimation
================

.. code:: python

   with gridpack.Session() as s:
       se = gridpack.StateEstimation(s, "input.xml")   # measurementList from the XML
       result = se.solve()
       print(result.converged)
       result.write()


Running in parallel
===================

Any of the above runs unchanged under MPI::

    mpiexec -np 4 python my_script.py

Keep a few buses of work per rank: a network too small for the rank count
hangs rather than failing.  The 9-bus case runs on up to three ranks; IEEE
14 and 118 are comfortable at four.  Print and write files from one rank,
but keep the collective calls listed above on every rank.


The command line
================

The same applications are available without writing Python::

    gridpack pf input_14.xml
    mpiexec -np 4 gridpack ca input_14_ca.xml --no-print-calcs

See :doc:`cli/index`.


Where to go next
================

``python/examples/`` in the source tree has five short scripts meant to be
read, one application each and then one study that chains three, with a
README saying what each shows and which rank counts were verified.  The
pages that follow are the reference for the low-level pybind11 objects.
