## Configuration Options

The contingency analysis application is configured via XML input files. Key options include:

### Contingency Specification

You can specify contingencies in two ways (or combine both):

**Option 1: Auto-generate N-1 contingencies**
```xml
<Contingency_analysis>
  <FullBranchN1>true</FullBranchN1>      <!-- N-1 for all branches -->
  <FullGeneratorN1>true</FullGeneratorN1> <!-- N-1 for all generators -->
</Contingency_analysis>
```

**Option 2: Use a contingency list file**
```xml
<Contingency_analysis>
  <contingencyList>contingencies.xml</contingencyList>
</Contingency_analysis>
```

**Option 3: Combine both** (auto-generate N-1 + custom N-K from file)
```xml
<Contingency_analysis>
  <FullBranchN1>true</FullBranchN1>
  <FullGeneratorN1>true</FullGeneratorN1>
  <contingencyList>custom_nk_contingencies.xml</contingencyList>
</Contingency_analysis>
```
When combined, duplicates from the file are automatically skipped.

### Other Options

| Option | Description | Default |
|--------|-------------|---------|
| `groupSize` | Number of MPI processes per contingency (parallelization) | 1 |
| `printCalcFiles` | Write detailed output for each contingency | true |
| `minVoltage` | Minimum voltage threshold for violations (p.u.) | 0.9 |
| `maxVoltage` | Maximum voltage threshold for violations (p.u.) | 1.1 |
| `qlim` | Enable reactive power limit enforcement (PV to PQ bus conversion) | false |
| `outputFormat` | `text` / `json` / `csv` / `csv_flat` / `csv_delta` | `text` |
| `outputFile` | Base name for output files | `ca_results` |
| `writeStats` | Emit StatBlock summary files (vmag.txt etc.). Set false to skip and avoid the per-case StatBlock work | true |
| `contingencyRating` | Loading% denominator across every CA output: `A`, `B`, or `C` with A→B→C fallback. Default `A` matches PW / PSS/E ACCC. `base_rate_mva` always uses rate-A | `A` |
| `monitorBranchesFile` | Path to a CSV allowlist (`from_bus,to_bus,ckt`). When set, overrides the area/kV gates (TARA / PSS/E convention) | (unset) |
| `monitorAreas` | Space-separated list of PSS/E area numbers. Branch is emitted if **either endpoint** is in the set | (unset) |
| `monitorKvMin` | Lower kV threshold; branch passes if `max(kv_from, kv_to) >= monitorKvMin` | 0 (unbounded) |
| `monitorKvMax` | Upper kV threshold; branch passes if `max(kv_from, kv_to) <= monitorKvMax` | 0 (unbounded) |

### Filtering csv_flat / csv_delta output

`csv_flat` and `csv_delta` emit per-(contingency, branch) rows. All filters
are optional — unset means "monitor everything". `monitorBranchesFile` is
authoritative when set; otherwise `monitorAreas` and the kV bounds AND
together.

```xml
<!-- (a) Monitor everything (default): no filter options set -->

<!-- (b) Curated allowlist -->
<monitorBranchesFile>monitor_branches.csv</monitorBranchesFile>

<!-- (c) Topology-based screening -->
<monitorAreas>11 12 19</monitorAreas>
<monitorKvMin>100.0</monitorKvMin>
<monitorKvMax>500.0</monitorKvMax>
```

`monitorBranchesFile` is a CSV of `from_bus,to_bus,ckt` rows (header
optional; `#` is a line comment). When set, area/kV options are ignored
and a warning is logged.

`monitorAreas` matches branches with either endpoint in the set (catches
tie-lines). `monitorKvMin/Max` gate on `max(kv_from, kv_to)` so a 138/13.8
step-down counts as 138.

`contingencyRating` (`A` | `B` | `C`, default `A`) sets the loading%
denominator for every CA output. Default `A` matches PowerWorld / PSS/E
ACCC. `base_rate_mva` always uses rate-A. A→B→C fallback if the picked
tier is zero/missing.

A complete annotated example is in
`src/applications/data_sets/input/ca/input_14_filters_example.xml` with a
sample monitor file `monitor_branches_14.csv` in the same directory.

### Output ordering

Rows are not sorted by contingency. The driver distributes contingencies
across MPI ranks and streams each rank's results to its own `.part` file;
rank 0 concatenates in rank order, so the final file is grouped by rank
and ordered by completion within each rank. Column 1 (`event_idx`)
preserves input-deck order — sort downstream if needed:

```bash
( head -1 my_run_delta.csv && tail -n +2 my_run_delta.csv | sort -t, -k1,1n ) > my_run_delta.sorted.csv
```


### Contingency File Format

See `contingencies_nk_example.xml` for examples of N-1, N-2, and N-3 contingency definitions.

The line-tag and generator-id elements accept PSS/E-aligned aliases for clarity:

| Element | Alias | Holds |
|---|---|---|
| `<contingencyLineNames>` | `<CKT>` | Branch circuit ID (PSS/E `CKT` field) |
| `<contingencyGenerators>` | `<GenID>` | Generator ID (PSS/E `ID` field) |

Either name works; mix-and-match within the same file is fine. Both legacy
files and new files using the PSS/E names continue to parse without changes.

---

## Advanced Features

### Automatic Slack Bus Transfer

When a generator contingency trips the slack bus generator, the application automatically transfers the slack bus role to the largest remaining online generator. This mimics commercial power flow tools (PSS/E, PowerWorld) behavior.

**Example output:**
```
Slack bus transferred from bus 69 to bus 80 (capacity: 577.0 MW)
```

After the contingency analysis completes, the slack bus is restored to its original location.

### Slack Capacity Check

After the power flow solves, the application checks if the slack bus generator output exceeds its Pmax rating. If the required generation exceeds capacity, the contingency is marked as failed with a warning message:

### Island Detection

The application detects network islands (disconnected portions) caused by branch contingencies:

- **Lone bus**: A bus with no active branches is marked as isolated
- **Island**: A group of buses disconnected from the main network

When isolation is detected, the isolated buses are excluded from the power flow solution, and a warning is added to the results. The main network portion is still solved if possible.

---

## Output Files

The contingency analysis calculation produces a number of files summarize the
results of the entire set of individual contingency simulations. Only results
for contingencies that ran to completion are included. Calculations that failed
either because of a numerical instability or because the calculations failed to
converge are not included in the results. The output files are described below.

### CSV outputs (`outputFormat=csv_flat` / `csv_delta`)

When `outputFormat` is set to `csv_flat` or `csv_delta`, the application writes
per-(contingency, branch) rows for downstream statistical analysis instead of
the aggregated `.txt` files described later in this section. All file names
below use the value of `outputFile` as a prefix; the default prefix is
`ca_results`.

**`<outputFile>_delta.csv`** *(`outputFormat=csv_delta` only)* — wide-form,
one row per (contingency, monitored branch) joining base + contingency state
on the same row. This is the format most downstream consumers prefer because
each row is self-contained (no separate base-case join needed).

| # | Column | Notes |
|---|---|---|
| 1 | `event_idx` | 0 = base case (only if a base row is emitted), 1..N = contingencies in input-deck order |
| 2 | `contingency` | contingency name from the input XML (`base_case` for the base) |
| 3 | `type` | `branch` or `generator` — what kind of contingency was tripped |
| 4–6 | `from_bus`, `to_bus`, `ckt` | Identity of the **monitored branch** in the row (not the tripped element) |
| 7–8 | `base_kv_from`, `base_kv_to` | Endpoint base kV |
| 9–10 | `area_from`, `area_to` | PSS/E area numbers |
| 11 | `base_rate_mva` | Always rate-A (PSS/E "normal" rating) |
| 12 | `cont_rate_mva` | Rating selected by `contingencyRating` (default A, with A→B→C fallback if zero/missing) |
| 13–14 | `base_p_mw`, `cont_p_mw` | Real-power flow before / after contingency |
| 15–16 | `base_q_mvar`, `cont_q_mvar` | Reactive-power flow before / after |
| 17–18 | `base_mva`, `cont_mva` | `sqrt(P² + Q²)` before / after |
| 19 | `base_loading_pct` | `base_mva / base_rate_mva × 100` |
| 20 | `cont_loading_pct` | `cont_mva / cont_rate_mva × 100` |
| 21–22 | `v_from_base`, `v_from_cont` | From-bus voltage magnitude (pu) before / after |
| 23–24 | `v_to_base`, `v_to_cont` | To-bus voltage magnitude (pu) before / after |
| 25–28 | `ang_from_base`, `ang_from_cont`, `ang_to_base`, `ang_to_cont` | Bus angles (deg) |
| 29–30 | `d_angle_base`, `d_angle_cont` | `ang_from − ang_to` before / after |
| 31 | `cont_event_facility` | Identifier of the **tripped element** in this contingency (e.g. `[area] from to ckt` for branch trips, `gen <bus> <id>` for gen trips) |

**`<outputFile>_flat.csv`** *(`outputFormat=csv_flat` only)* — long-form,
one row per (case, branch). Columns: `event_idx, contingency, from_bus,
to_bus, ckt, p_from_mw, q_from_mvar, mva_from, rate_mva, loading_percent,
viol, v_from_pu, v_to_pu, ang_from_deg, ang_to_deg`. `rate_mva` is rate-A
on `event_idx=0` rows, the configured `contingencyRating` on contingency
rows.

**`<outputFile>_buses.csv`** *(both `csv_flat` and `csv_delta`)* — bus
metadata sidecar so the per-branch files can stay narrow. Columns:
`bus_id, bus_name, base_kv, area, zone, owner, area_name, zone_name,
owner_name`.

**`<outputFile>_convergence.csv`** *(every `outputFormat`)* — one row per
contingency. Columns: `event_idx, contingency, type, converged,
iterations, final_tolerance, max_p_bus, max_p_mismatch, max_q_bus,
max_q_mismatch, status_code`. `converged` is `true` iff `status_code ==
"OK"`; `status_code` is one of `OK` / `ISLANDED` / `NO_SLACK` /
`DIVERGED` / `SLACK_OVERLOAD`. Failed rows appear here even though
they're omitted from `_delta.csv` / `_flat.csv`.

When monitor filters are active, the data-row count of `_delta.csv` /
`_flat.csv` equals `|monitored branches| × |converged contingencies|`.

---

### Aggregated `.txt` outputs (`writeStats=true`, default)

The remaining files in this section are produced by the StatBlock summary
pipeline, controlled by the `writeStats` option (default `true`). They
contain per-element statistics (mean / RMS / min / max) aggregated **across
all contingencies**, not per-contingency rows. Set `writeStats=false` to
skip them when csv_flat / csv_delta output is sufficient.

The set of files emitted is fixed; their names are not configurable. With
`writeStats=true` you get: `vmag.txt`, `vmag_mm.txt`, `vang.txt`,
`vang_mm.txt`, `pgen.txt`, `pgen_mm.txt`, `qgen.txt`, `qgen_mm.txt`,
`pflow.txt`, `pflow_mm.txt`, `qflow.txt`, `qflow_mm.txt`, `perf_mm.txt`,
`perf_sum.txt`, `line_flt_cnt.txt`. The file `pq_change_cnt.txt` is also
written when `qlim=true`. Per-contingency convergence/status is reported
via the `_convergence.csv` sidecar described above, which is written for
every `outputFormat`.

**vmag.txt**: This file contains the average value of the voltage magnitude for
non-PV buses. It also contains the RMS fluctuations of the voltage magnitude
with respect to the voltage average and also with respect to the base case. The
columns in this file are

column 1: row index

column 2: bus ID

column 3: average voltage magnitude over all contingencies

column 4: RMS voltage magnitude fluctuations with respect to the average value

column 5: RMS voltage magnitude fluctuations with respect to the base case value

**vmag\_mm.txt**: This file contains the minimum and maximum values of the
voltage magnitude over all contingencies. It also contains the contingency index
for which the minimum or maximum value occurs.

column 1: row index

column 2: bus ID

column 3: base case value

column 4: minimum value over all contingencies

column 5: maximum value over all contingencies

column 6: deviation of minimum value from base case

column 7: deviation of maximum value from base case

column 8: contingency index of minimum value

column 9: contingency index of maximum value

**vang.txt**: This file has the same structure as vmag.txt, except that the
stored values all represent the phase angle at each bus. PV buses are included
in this data.

**vang\_mm.txt**: This file has the same structure as vmag\_mm.txt, except that the
stored values all represent the phase angle at each bus. PV buses are included
in this data.

**pq\_changed\_cnt.txt** This file is only created if the qlim flag is
set to true in the input file. It counts the number of times a PV bus is
changed to a PQ bus during the simulation.

column 1: row index

column 2: bus ID

column 3: number of contingencies where PV bus changed to PQ bus

**pgen.txt**: This file contains the average value of the real power for each
generator in the system. It also contains the RMS deviations of the real power
fluctuations with respect to the average and also with respect to the base case.

column 1: row index

column 2: bus ID

column 3: 2 character generator ID

column 4: average value of real power for each generator

column 5: RMS fluctuation in the real power with respect to average value

column 6: RMS fluctuation in the real power with respect to base case value

**pgen\_mm.txt**: This file contains the minimum and maximum values of the real
power for each generator, as well as the contingency index for which those
values occured.

column 1: row index

column 2: bus ID

column 3: 2 character generator ID

column 4: real power generation for the base case

column 5: minimum value of real power generation

column 6: maximum value of real power generation

column 7: deviation of minimum value from base case

column 8: deviation of maximum value from base case

column 9: contingency index of minimum value

column 10: contingency index of maximum value

**qgen.txt**: This file has the same structure as pgen.txt, except that the
stored values represent the reactive power at each generator.

**qgen\_mm.txt**: This file has the same structure as pgen\_mm.txt, except that
the stored values represent the reactive power at each generator.

**pflow.txt**: This file contains the average value of the real power flow at
the "from" bus for each line element. It also contains the RMS fluctuations with
respect to the average value and the base case.

column 1: row index

column 2: bus ID for "from" bus

column 3: bus ID for "to" bus

column 4: 2 character line ID

column 5: average value of real power flow for "from" bus

column 6: RMS fluctuations of real power flow with respect to average value

column 7: RMS fluctuations of real power flow with respect to base case value

**pflow\_mm.txt**: This file contains the minimum and maximum values of the real
power flow at the "from" bus for each line element. It also contains the
contingency indices at which these values occur.

column 1: row index

column 2: bus ID for "from" bus

column 3: bus ID for "to" bus

column 4: 2 character line ID

column 5: base case value of real power flow

column 6: minimum value of real power flow

column 7: maximum value of real power flow

column 8: deviation of minimum value from the base case value

column 9: deviation of maximum value from the base case value

column 10: minimum allowable value of power flow (-rate A parameter)

column 11: maximum allowable value of power flow (rate A parameter)

column 12: contingency index of minimum value

column 13: contingency index of maximum value

**qflow.txt**: This file has the same structure as pflow.txt except that all
values are for the reactive power flow at the "from" bus.

**qflow_mm.txt**: This file has the same structure as pflow\_mm.txt except that
all values are for the reactive power flow at the "from" bus.

**perf\_mm.txt**: This file contains minimum and maximum values of the real
power flow performance index (this is the absolute value of the complex power on
each line divided by the line rating value, squared). It also contains the
contingency indices for which these values occur.

column 1: row index

column 2: bus ID for "from" bus

column 3: bus ID for "to" bus

column 4: 2 character line ID

column 5: base case value of the performance index

column 6: minimum value of the performance index

column 7: maximum value of the performance index

column 8: deviation of minimum value from the base case

column 9: deviation of maximum value from the base case

column 10: contingency index at which the minimum value occurs

column 11: contingency index at which the maximum value occurs

**perf\_sum.txt**: This file contains the sum of the performance index over all
lines for each contingency. It also contains the average value of the
performance index over all lines for each contingency.

column 1: contingency index (0 is the base case)

column 2: sum of performance index over all lines

column 3: average value of performance index over all lines

**line\_flt\_cnt.txt**: This file contains the total number of faults found on
each line for all contingencies.

column 1: row index

column 2: bus ID for "from" bus

column 3: bus ID for "to" bus

column 4: 2 character line ID

column 5: total number of contingencies that result in a fault on this line
