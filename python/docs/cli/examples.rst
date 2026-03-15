Usage Examples
==============

This chapter presents representative usage scenarios for each
subcommand, illustrating common workflows and option combinations.

Power Flow
----------

Solve the IEEE 14-bus system using the default Newton-Raphson solver::

    gridpack powerflow IEEE_14_bus.xml

Use the PETSc nonlinear solver and write results to a file::

    gridpack pf IEEE_14_bus.xml --solver nl --output pf_results.txt

Export the solved case to PSS/E v33 format::

    gridpack pf IEEE_14_bus.xml --export-psse v33 IEEE14_solved.raw

Run on four MPI processes for a larger network::

    mpiexec -np 4 gridpack pf large_network.xml


Dynamic Simulation
------------------

Perform a transient stability study on the IEEE 39-bus system::

    gridpack dsf IEEE39bus_ds.xml

Run in parallel with output redirected to a file::

    mpiexec -np 8 gridpack dsf IEEE39bus_ds.xml --output dsf_results.txt

Suppress console output during a batch run::

    gridpack dsf IEEE39bus_ds.xml --quiet --no-timer


State Estimation
----------------

Run state estimation using measurements defined in the XML configuration::

    gridpack se IEEE14_se.xml

The measurement data file is specified by the ``measurementList``
parameter within the ``<State_estimation>`` block of the XML
configuration. Results can be redirected to a file::

    gridpack se IEEE14_se.xml --output se_results.txt


Contingency Analysis
--------------------

Analyze all contingencies defined in the contingency list file::

    gridpack ca IEEE_14_bus_ca.xml

Tighten voltage limits and suppress per-contingency output files::

    gridpack ca IEEE_14_bus_ca.xml --vlimits 0.95 1.05 --no-print-calcs

Run contingencies in parallel across four MPI processes::

    mpiexec -np 4 gridpack ca large_network_ca.xml


HADREC
------

Run a HADREC simulation with default settings::

    gridpack hadrec IEEE39bus_hadrec.xml

Suppress output for automated testing::

    gridpack hadrec IEEE39bus_hadrec.xml --quiet


EMT Simulation
--------------

Run an EMT simulation::

    gridpack emt emt_config.xml

Run in parallel::

    mpiexec -np 4 gridpack emt emt_config.xml


Workflow Integration
--------------------

The CLI integrates naturally with shell scripts and automation
pipelines. The following example runs a base-case power flow followed
by contingency analysis on the same network:

.. code-block:: bash

    #!/bin/bash
    set -e

    echo "Solving base-case power flow..."
    gridpack pf network.xml --output base_case.txt

    echo "Running contingency analysis..."
    mpiexec -np 4 gridpack ca network_ca.xml --no-print-calcs

    echo "Analysis complete."

Output Redirection
~~~~~~~~~~~~~~~~~~

The ``--output`` flag redirects structured application output (bus
voltages, branch flows, etc.) to a file. Diagnostic messages and
solver progress information are still printed to the console. To
capture all output, use shell redirection::

    gridpack pf input.xml > all_output.txt 2>&1
