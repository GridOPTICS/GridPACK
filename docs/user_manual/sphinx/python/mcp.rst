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

The server needs Python 3.10 or later and the ``mcp`` extra; the ``plot``
extra adds matplotlib for ``plot_job``::

    pip install 'gridpack[mcp,plot]'
    gridpack-mcp doctor

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
``gridpack-mcp config`` prints both registrations, the ``claude mcp add``
line and the JSON block, for the options you give it, with a ``PATH`` that
covers the server, ``gridpack`` and ``mpiexec``::

    gridpack-mcp config --parallel 3 --cases ~/gridpack-cases

Options
~~~~~~~

``--timeout SECONDS``
    Stop any run after this long; the default is 600.  Clients cannot
    change it.

``--parallel N``
    Run up to ``N`` jobs at a time (default 1).  Studies fan out through
    this; a job's ranks times ``N`` may not exceed the machine's cpus.

``--cases DIR``
    Also offer the input files under ``DIR``, and let ``create_case`` build
    new ones from RAW and DYR files there.  May be repeated.

``--output DIR``
    Where ``save_figure``, ``save_report`` and ``export_data`` write.
    Without it they write under the server's temporary directory, which is
    gone when the server exits, and their replies say so.

``gridpack-mcp doctor`` checks all of this (Python, ``mcp``, the
``gridpack`` CLI, ``mpiexec``, matplotlib, the cases it would offer) and
ends with ``ready`` or ``not ready``.

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

* **Cases from a RAW file** are made by ``create_case``: give it a RAW,
  and optionally a DYR, under a ``--cases`` directory and it writes an
  input XML beside them with default power flow and dynamic simulation
  settings, every generator watched and observed and every bus observed,
  and no events.  Without a DYR of your own, ``create_dyr`` writes one with
  generic models for every generator in the RAW, ``classical`` (GENCLS) or
  ``detailed`` (GENROU with a SEXS exciter and a TGOV1 governor), at
  textbook parameters: fine for teaching, screening and comparing
  scenarios, not a model of any real machine.
* **State estimation cases** come from ``create_measurements``: it takes a
  finished power flow job as the truth and writes a measurement file with
  chosen types, noise and planted bad data (see :ref:`mcp-measurements`).
* **Derived cases** are made by ``modify_case`` (below) and live under the
  server's work directory as ``derived/<name>/<name>.xml`` until the
  server exits.  An ``edits.json`` beside each records how it was made.

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

   * - ``job_series``
     - Chosen channels of a dynamic run, by name or glob, decimated to a
       point cap.
   * - ``job_report``
     - A short text verdict: convergence, totals, voltage range, violations
       and overloads for a power flow; events, speed deviation, damping,
       rotor angle spread and a stability verdict for a dynamic run.
   * - ``plot_job``
     - A PNG: generator speeds or any channels with fault windows shaded
       and trips marked, or a voltage or angle profile.  Needs the
       ``plot`` extra.
   * - ``compare_jobs``
     - What changed between two runs of the same kind.
   * - ``create_case``, ``create_dyr``, ``modify_case``,
       ``create_measurements``
     - New cases from RAW/DYR files, generic DYR models, derived cases with
       edits, events, observations and model parameters, and state
       estimation cases from a power flow (next sections).
   * - ``run_study``, ``study_status``, ``study_result``
     - Load sweeps and dynamic N-1 screens, one job per point, ranked.
   * - ``save_figure``, ``save_report``, ``export_data``
     - Files for papers and readers: vector figures with a recipe,
       self-contained HTML or Markdown reports, CSV or JSON data (see
       :ref:`mcp-files`).

Jobs are kept in memory and run ``--parallel`` at a time; the server stops
any that are still running when it exits.  EMT runs return status only.

Deriving cases
--------------

``modify_case(case, edits, name, replace_events)`` writes a copy of the
case with its RAW records edited and its events changed, and lists it.
Only the named fields of the named records change; the rest of the RAW is
kept byte for byte, for RAW versions 23 to 36.  Every edit is checked
before anything is written, and the reply echoes each one with the value
before and after.

Network edits, in engineering units (ids default to the only one at the
bus; circuit ids must match the RAW):

* ``set_load_active_power {bus, mw}``, ``set_load_reactive_power {bus,
  mvar}``, ``scale_load {bus, factor}``, ``set_load_in_service {bus,
  in_service}``, ``scale_all_loads {factor}``
* ``set_generator_active_power {bus, mw}``,
  ``set_generator_voltage_setpoint {bus, vm_pu}``,
  ``set_generator_in_service {bus, in_service}``,
  ``scale_all_generation {factor}``
* ``set_shunt {bus, mvar, mw}``, ``set_shunt_in_service {bus, in_service}``
* ``set_branch_in_service {from, to, ckt, in_service}``
* ``set_param {kind, bus | from, to, ckt, id, field, value}`` with a PSS/E
  RAW column name, for anything else.

Event edits, times in seconds on the case's time step:

* ``add_fault {from, to, start, duration}``: a solid fault at the from bus,
  cleared after ``duration`` (GridPACK's ``faultEvent``).
* ``add_bus_fault {bus, start, duration, g_pu, b_pu}``.
* ``trip_branch {from, to, ckt, time, in_service}`` and
  ``trip_generator {bus, id, time, in_service}``.
* ``add_observations {generators, buses, bus_frequency}`` with ``"all"`` or
  bus lists: what HADREC records.
* ``set_dyr_param {bus, id, model, field, value}``: a dynamic model
  parameter in the case's DYR, with the field names GridPACK's parsers use
  (GENCLS H, D, XDP; GENROU TDOP, TDOPP, TQOP, TQOPP, H, D, XD, XQ, XDP,
  XQP, XDPP, XL, S1, S12; SEXS TA_TB, TB, K, TE, EMIN, EMAX; TGOV1 R, T1,
  VMAX, VMIN, T2, T3, DT) or ``p3``, ``p4``, ... by position for any other
  model.

Events are appended to the base case's own unless ``replace_events`` is
true.  HADREC needs a fault as the first event and refuses a case with no
observations.  Not editable yet: transformer taps and status in v33 and
later files, and switched shunts.

.. _mcp-measurements:

Measurements for state estimation
---------------------------------

``create_measurements(job_id, ...)`` takes a finished, converged power flow
job as the truth and writes a state estimation case: a measurement file
with the chosen types (VM, VA in degrees, PI, QI, PIJ, QIJ, PJI, QJI in
p.u. on the case base; by default every bus voltage, injection and
from-end flow), Gaussian noise per type whose standard deviation is also
each meter's ``Deviation``, ``coverage`` to keep a fraction of the meters,
and planted gross errors by meter or ``random_bad`` of them at
``bad_sigma`` standard deviations, all repeatable with ``seed``.  Running
``run_state_estimation`` on the case and then ``job_report`` adds the
estimate's RMS and largest errors against the truth and which planted bad
data the detector caught, missed, or raised false alarms on.
``noise_scale`` 0 writes exact values but keeps the nominal deviations, as
the estimator weights meters by them.

Studies
-------

``run_study(case, kind, ...)`` writes one derived case per point, runs one
job each, ``--parallel`` at a time, and returns a study id for
``study_status`` and ``study_result``.

``load_sweep``
    ``factors`` such as ``[1.0, 1.1, 1.2]`` scale every load (and every
    generator with ``scale_generation``) and run power flows.  The table
    gives load, generation and losses in MW, the lowest voltage and its
    bus, violations and peak loading per factor.

``dynamic_n1``
    Every branch and transformer, or the ``branches`` given as
    ``"from-to"`` or ``"from-to ckt"``, gets a fault at its from bus at
    ``fault_start`` cleared after ``fault_duration`` by tripping it, with
    the base case's events replaced and every machine watched, and a
    dynamic simulation to ``end_time``.  The table ranks outages most
    severe first: ``unstable`` means the rotor angle spread between
    machines, each measured from its own start, passed 180 degrees, or a
    machine was still more than 2% off speed at the end.  Treat the
    180-degree rule as a first-swing screen; a large swing that recovers
    can trip it, so read the final spread and speed columns too.

.. _mcp-files:

Files for papers and readers
----------------------------

``plot_job`` answers in the chat.  Three tools write files instead, into
the ``--output`` directory:

``save_figure(job_ids, ...)``
    A figure of one or more finished jobs as ``pdf`` (default), ``svg`` or
    ``png``, with a JSON recipe beside it naming the jobs, their cases and
    every setting, so the figure can be regenerated after a run is redone.
    ``style`` ``paper`` drops the title (the caption belongs in the paper),
    uses print sizes with embedded fonts, labels axes with units, and gives
    every series a dash pattern as well as a color so it survives
    grayscale; ``screen`` is the chat style.  ``width`` is ``single`` (3.5
    in), ``double`` (7.16 in) or a number of inches; ``height`` is optional.
    Several jobs overlay: one color per job, named by ``labels``, and a dash
    per channel, or per job when they share one channel.  ``kind`` and
    ``channels`` are as for ``plot_job``; ``font`` is ``serif`` or ``sans``.

``save_report(job_id | study_id, ...)``
    One file a reader can open without the server: the case and, for a
    derived case, the edits that made it; the events; the solver settings;
    the ``job_report`` verdict; the result tables; and figures, inline in
    ``html`` or as PNG files beside a ``md`` report.  A study report leads
    with the ranked table and shows the three worst outages.

``export_data(job_id, ...)``
    The time series of a dynamic run at full resolution, with chosen
    channels, or the bus, branch, generator and contingency tables of a
    static run, as ``csv`` or ``json``, for your own plotting.

Names are given with ``name`` or taken from the job; a file is never
overwritten unless ``overwrite`` is true.

Prompts and study workflow
--------------------------

The server offers three prompt templates, ``stability_study``,
``n1_screen`` and ``voltage_study``, that lead a client through a study in
the right order.  The same workflow, for any client:

* Start with ``list_cases`` and ``describe_case``: note ``max_ranks``, the
  solver settings (time step and end time) and what a dynamic case watches
  and observes.
* Event times must be multiples of the case's time step and inside its
  length.  When the base case already has a fault, pass ``replace_events``
  so only the new events run.  To clear a fault by opening the line, add a
  ``trip_branch`` at ``start + duration``.
* If a run tool returns ``running``, poll ``job_status``; do not resubmit.
  ``channels`` 0 on a dynamic job means the case watches nothing.
* Read ``job_report`` first, then ``plot_job`` for figures in the chat,
  ``job_series`` with a channel list and ``max_points`` to quote values, and
  ``compare_jobs`` for a base against a derived run.
* For a paper, ``save_figure`` with the ``paper`` style and the column
  ``width``, captioned from ``job_report``'s numbers; start the server with
  ``--output DIR`` so the files outlive it.
* In a ``dynamic_n1`` table, outages that island a generator (its step-up
  branch) rank worst, and rightly so.
* Reports judge voltages against 0.9 to 1.1 p.u. and overloads against
  rate A; a case without ratings cannot be judged for overloads.

.. autoclass:: gridpack.mcp.tools.Tools
   :members: list_cases, describe_case, create_case, create_dyr,
             modify_case, create_measurements, run_powerflow, run_contingency, run_dynamic_sim,
             run_state_estimation, run_hadrec, run_emt, job_status,
             job_result, job_log, job_series, job_report, plot_job,
             compare_jobs, run_study, study_status, study_result,
             save_figure, save_report, export_data

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

    **You:** Screen every branch outage of the 9-bus case over one second
    and show me the worst.

    ``run_study(case="input/ds/input_9b3g.xml", kind="dynamic_n1",
    end_time=1.0)`` then ``study_result(study_id="s33757c6e")`` ::

        {"finished": 9, "pending": 0, "failed": [], "table": [
          {"branch": "3-6 BL", "unstable": true, "angle_spread_peak_deg": 165.4,
           "peak_speed_deviation": 0.0409, "final_speed_deviation": 0.0409, ...},
          {"branch": "1-4 BL", "unstable": true, ...},
          {"branch": "2-8 BL", "unstable": true, ...},
          {"branch": "8-9 BL", "unstable": false, ...}, ...]}

    ``job_report(job_id=...)`` for the first row ::

        Dynamic simulation of 0.99 s at 0.01 s steps, 100 rows, 24 channels.
        Event: fault at bus 3 from 0.1 to 0.15 s (50 ms).
        Event: branch out at 0.15 s (line 3 6).
        Peak speed deviation 4.092% on generator 1 at bus 3; final 4.092%.
        Verdict: large speed excursion; speed still 4.1% off at the end:
        machine running away.

    Opening 3-6, 1-4 or 2-8 islands a generator, so those three rank worst.
    ``plot_job(job_id=...)`` returns the speed curves as an image.

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
* **HADREC finished but recorded only time**: the case has no
  observations; derive one with ``add_observations``.  ``run_hadrec`` now
  refuses such a case with that advice.
* **plot_job says to install matplotlib**: ``pip install 'gridpack[plot]'``
  in the server's environment.
* **run_study refuses for cpus**: lower ``ranks`` or start the server with
  a smaller ``--parallel``.
