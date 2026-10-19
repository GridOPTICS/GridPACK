---
name: gridpack-study
description: Run power system studies through the GridPACK MCP server (gridpack-mcp) - pick or create a case, change loads, generators, shunts and branches, add faults and trips, run power flow, dynamic simulation, HADREC, contingency screens, N-1 studies and load sweeps, and read back reports, plots and comparisons. Use whenever the user wants to simulate, screen, stress or compare a power system case with GridPACK.
---

# GridPACK studies over MCP

The gridpack server runs GridPACK in subprocesses and returns data. These
notes say how to use its tools well; the tool descriptions say what each
argument means.

## Start every study the same way

1. `list_cases` (optionally with `app`) to find the case, then
   `describe_case`. Note three things: `max_ranks`, the solver `settings`
   (time step and end time for dynamics), and `records` for dynamic cases
   (how many generators are watched, what HADREC observes).
2. Never ask for more ranks than `max_ranks`. A 9-bus case runs on at most
   3 ranks; split further, a rank gets no buses and every rank spins. Say
   this plainly if the user asks for more.
3. The user's own RAW and DYR files must sit under a directory the server
   was started with (`--cases DIR`); `create_case` builds the input from
   them. Cases are named `<DIR name>/<path>`. With no DYR, `create_dyr`
   writes generic models (`classical` or `detailed`) for every generator;
   always say these are textbook parameters, not the real machines.
4. For state estimation, `create_measurements` on a finished power flow
   job writes meters with noise and planted bad data; `job_report` on the
   estimate then compares it with the truth.

## Changing a case

`modify_case` derives a new case; the base is never touched. Edits are a
list of ops in engineering units: `set_load_active_power {bus, mw}`,
`scale_load {bus, factor}`, `scale_all_loads {factor}`,
`set_generator_active_power {bus, mw}`, `set_generator_in_service`,
`set_shunt {bus, mvar}`, `set_branch_in_service {from, to, ckt, in_service}`
and more. Ids default to the only one at the bus; circuit ids must match the
RAW exactly, and the error message lists the choices when they do not.

- The reply echoes each edit with `was` and `now`; quote those numbers.
- Derived cases live only while the server runs. `edits.json` beside each
  one records how it was made, so a study is reproducible from the base.
- `set_dyr_param {bus, model, field, value}` tunes a dynamic model
  parameter (GENROU `H`, `XDP`; SEXS `K`; TGOV1 `R`; or `p3`, `p4`, ... by
  position). Transformer taps cannot be edited yet; say so.

## Faults and trips

Events go in the same `modify_case` call: `add_fault {from, to, start,
duration}` (a solid fault at the from bus, cleared after `duration`),
`add_bus_fault`, `trip_branch {from, to, time}`, `trip_generator {bus, time}`.

- Every time must be a multiple of the case's time step, inside its length.
- Events are appended to the base case's own events. When the base already
  has a fault, pass `replace_events: true` so only the user's event runs.
- To clear a fault by opening the line, add the fault and a `trip_branch`
  at `start + duration`.
- HADREC initialises from the first event, which must be a fault, and
  records only the case's observations. Most bundled cases have none: add
  `{"op": "add_observations", "generators": "all", "buses": "all"}` or
  `run_hadrec` will refuse with that advice.

## Running and waiting

The run tools wait up to 10 s and return the job status. If `state` is
still `running`, poll `job_status`; do not resubmit. `channels` on a
dynamic job says how many quantities were recorded; 0 means the case
watches nothing. A `failed` job's reason is in `job_status` and the tail of
`job_log` with `stream: "stderr"`.

## Reading results

1. `job_report` first. It is written for a reader: convergence, totals,
   voltage range, violations and overloads for power flow; events, peak and
   final speed deviation, damping ratio, rotor angle spread and a verdict
   for dynamics. Pass its verdict on, with the numbers.
2. `plot_job` for the standard figures: generator speeds (default), any
   channels by glob such as `*_PowerAngle`, or a voltage profile for power
   flow. It returns an image; fault windows are shaded and trips marked.
3. `job_series` with a channel list and `max_points` when you want to draw
   your own chart or quote values; never pull a whole table.
4. `compare_jobs` base against derived: it names the bus with the largest
   voltage change, lists branches present in one run only, and for dynamics
   the first time the runs differ and the channel that moved most.

## Figures and reports for a paper

`plot_job` is for the chat. When the user wants a figure to keep, use
`save_figure` with the paper style (the default): PDF or SVG, `width`
`single` or `double` for the journal column, `labels` naming the runs when
several overlay, and no title, because the caption goes in the paper. Say
where the file and its JSON recipe landed; if the reply carries a note
that the output is temporary, tell the user to start the server with
`--output DIR`. Write the caption yourself from `job_report`'s numbers, and
offer the LaTeX figure block:

```latex
\begin{figure}[t]\centering
\includegraphics[width=\columnwidth]{fig1.pdf}
\caption{...}\label{fig:...}\end{figure}
```

`save_report` writes a self-contained HTML (or Markdown) report of a job or
a study for a reader without the server; `export_data` writes the plotted
series or the result tables as CSV for the user's own tooling.

## Studies

`run_study` fans out one derived case and one job per point and ranks them:

- `kind: "load_sweep"` with `factors` scales every load (and generation with
  `scale_generation`) and runs power flows; the table gives load, losses,
  lowest voltage and its bus, violations and peak loading per factor.
- `kind: "dynamic_n1"` faults and trips every branch (or a `branches`
  list) and runs dynamics; the table ranks outages most severe first.

Poll `study_status` until `pending` is 0, then `study_result`. Explain the
ranking honestly: `unstable` means the rotor angle spread between machines
(each from its own start) passed 180 degrees, or a machine was still more
than 2% off speed at the end. The 180-degree rule is a first-swing screen;
a large swing that recovers can trip it, so look at `angle_spread_final_deg`
and the speed columns before calling an outage lost. Outages that island a
generator (its step-up branch) rank worst and are real.

## Limits to state up front

- EMT runs return status only, no results.
- `max_ranks` times `--parallel` jobs is capped at the machine's cpus;
  `run_study` says so when the budget is exceeded.
- Jobs are stopped after the server's `--timeout` (default 600 s).
- Reports judge voltages against 0.9 to 1.1 p.u. and overloads against
  rate A; cases without ratings cannot be judged for overloads.
