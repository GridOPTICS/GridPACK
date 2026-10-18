MCP Server
==========

``gridpack-mcp`` lets an AI assistant, or any other `Model Context
Protocol <https://modelcontextprotocol.io>`_ client, run GridPACK studies:
list the available cases, start a power flow, contingency screen, dynamic
simulation, state estimation, HADREC or EMT run, and read the results back
as structured data.

Each run is a ``gridpack`` command (:doc:`cli/index`) in its own
subprocess and working directory, so a crash in GridPACK never takes the
server down, and parallel runs go through ``mpiexec`` as usual.  The
server talks to its client over stdio.

Installation
------------

The server needs Python 3.10 or later and the ``mcp`` extra::

    pip install 'gridpack[mcp]'

``gridpack`` and, for runs on more than one rank, ``mpiexec`` must be on
the server's ``PATH``.  The bundled cases are read from the
``src/applications/data_sets`` directory of a GridPACK source tree, found
through ``GRIDPACK_TEST_DATA``, then ``GRIDPACK_ROOT_DIR``, then the tree
the package was installed from.  Without one, start the server with
``--cases`` (see :ref:`mcp-own-cases`).

The image built from ``Dockerfile.pygridpack`` (:doc:`cli/docker`) already
contains all of this; nothing else is needed on the host but Docker.

Starting the server
-------------------

The client starts the server.  With Claude Code, register it once::

    claude mcp add gridpack -- gridpack-mcp

If ``gridpack`` and ``mpiexec`` are not on the ``PATH`` Claude Code passes
on, set it for the server, for example for a virtual environment and
Homebrew::

    claude mcp add gridpack \
        -e PATH="$HOME/.venvs/pygridpack/bin:/opt/homebrew/bin:$PATH" \
        -- ~/.venvs/pygridpack/bin/gridpack-mcp

From the Docker image::

    claude mcp add gridpack -- docker run -i --rm pygridpack gridpack-mcp

Other clients take the same command in their configuration file, usually
in this form:

.. code-block:: json

    {"mcpServers": {"gridpack": {
        "command": "docker",
        "args": ["run", "-i", "--rm", "pygridpack", "gridpack-mcp"]}}}

``claude mcp list`` shows whether the server starts and connects.

Options
~~~~~~~

``--timeout SECONDS``
    Stop any run after this long; the default is 600.  Clients cannot
    change it.

``--cases DIR``
    Also offer the input files under ``DIR``.  May be repeated.

.. _mcp-own-cases:

Cases
-----

Runs take a case name from ``list_cases`` and nothing else, so a client
cannot point the server at arbitrary files.

* **Bundled cases** are the inputs under ``data_sets/input`` and
  ``data_sets/emt``, named by their path there, such as
  ``input/powerflow/input_14.xml``.  Those whose data files are missing,
  or that no application can run, are left out.
* **Your own cases** come from the directories given with ``--cases``.
  Each input XML under ``DIR`` is named ``<DIR name>/<path>``::

      gridpack-mcp --cases ~/models        # models/feeder.xml, models/sub/dyn.xml

  In Docker, mount the directory; it can be read-only, since runs work in
  their own directories::

      docker run -i --rm -v ~/models:/models:ro pygridpack gridpack-mcp --cases /models

  The RAW, DYR and other files an XML names are looked for beside it, then
  in the server's working directory, then in ``DIR/raw``, ``DIR/dyr`` and
  the other ``data_sets`` subdirectories.  Contingency and measurement
  lists in the directory are not offered as cases; an input that cannot be
  used is reported on the server's stderr with the reason.  An EMT case
  uses a ``.petscrc`` beside it if there is one, else the bundled one.

Each case lists the applications it can run, judged from its blocks, its
bus count and the most ranks it may use.  A listed case can still fail in
GridPACK; ``job_status`` then says why.  A case with fewer than 9 buses
runs serially, a 9-bus case on at most three ranks (split further, a rank
gets no buses and every rank spins), and larger cases on up to as many
ranks as the machine has cores.  Requests for more are refused before
anything starts.

Tools
-----

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Tool
     - Purpose
   * - ``list_cases``
     - Available cases, optionally for one application.
   * - ``describe_case``
     - A case's applications, data files and current solver settings.
   * - ``run_powerflow``, ``run_contingency``, ``run_dynamic_sim``,
       ``run_state_estimation``, ``run_hadrec``, ``run_emt``
     - Start a run with optional overrides (tolerance and iterations, or
       time step and end time) and a rank count; wait up to 10 seconds,
       then return the job's status.
   * - ``job_status``
     - queued, running, done or failed, with the exit code and any error.
       A run that did not converge is done, with ``converged`` false.
   * - ``job_result``
     - The full results of a finished job: what ``gridpack <command>
       --json`` writes (:doc:`cli/commands`).
   * - ``job_log``
     - The last lines of a job's stdout or stderr.

Jobs are kept in memory and run one at a time; the server stops any that
are still running when it exits.

.. autoclass:: gridpack.mcp.tools.Tools
   :members: list_cases, describe_case, run_powerflow, run_contingency,
             run_dynamic_sim, run_state_estimation, run_hadrec, run_emt,
             job_status, job_result, job_log

Example
-------

A session in Claude Code, with the tool calls it makes and what they
return (abridged):

    **You:** Which cases can run state estimation?

    ``list_cases(app="se")`` ::

        [{"case": "input/se/input_118.xml", "apps": ["se"],
          "network": "IEEE118.raw", "buses": 118, "max_ranks": 12},
         {"case": "input/se/input_14.xml", "apps": ["se"],
          "network": "IEEE14.raw", "buses": 14, "max_ranks": 12}]

    **You:** Run a power flow on input/powerflow/input_14.xml and show bus 5.

    ``run_powerflow(case="input/powerflow/input_14.xml")`` ::

        {"id": "0c65bb165cd5", "app": "powerflow", "ranks": 1,
         "state": "done", "exit_code": 0, "error": null, "converged": true, ...}

    ``job_result(job_id="0c65bb165cd5")``, bus 5 of ``buses`` ::

        {"busId": 5, "name": "BUS-5", "voltage": 1.0195138598224722,
         "angle": -8.773853898144809, "pLoad": 7.6, "qLoad": 1.6, ...}

    **You:** Now a dynamic simulation of input/ds/input_9b3g.xml on 4 ranks.

    ``run_dynamic_sim(case="input/ds/input_9b3g.xml", ranks=4)`` ::

        Error executing tool run_dynamic_sim: input/ds/input_9b3g.xml
        (9 buses) runs on at most 3 ranks

Troubleshooting
---------------

* **A run failed**: ``job_status`` gives the reason; ``job_log`` with
  ``stream="stderr"`` shows what GridPACK printed.
* **"could not start"**: ``gridpack`` or ``mpiexec`` is not on the
  server's ``PATH``; set it when registering the server.
* **The server exits at once with "no cases"**: it found no ``data_sets``
  directory and was given no ``--cases``; set ``GRIDPACK_TEST_DATA`` to a
  GridPACK source tree.
* **An EMT run aborts at start**: it needs a ``.petscrc`` in its working
  directory; the server supplies one, so check a ``--cases`` directory's
  own copy.
