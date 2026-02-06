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

### Contingency File Format

See `contingencies_nk_example.xml` for examples of N-1, N-2, and N-3 contingency definitions.

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

```
WARNING: Slack bus 80 generator output (475.3 MW) exceeds capacity (400.0 MW)
Insufficient generation capacity for contingency GN_69_1
```

This ensures realistic results - a contingency that requires more generation than available capacity is properly flagged as a failure.

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

**success.txt**: This file summarizes the results of each contingency and
reports 1) whether the contingency calculation successfully ran to completion,
2) whether a violation was found (bus, branch, or both), and 3) whether any
buses were isolated.

Example output:
```
contingency: 1 success: true violation: none
contingency: 2 success: true violation: branch
contingency: 3 success: true violation: none warning: isolated
contingency: 4 success: false
```

- `success: true` - Power flow converged and slack capacity is within limits
- `success: false` - Power flow failed, island detected, or slack capacity exceeded
- `violation: none/bus/branch` - Whether voltage or thermal limits were violated
- `warning: isolated` - One or more buses were isolated (lone bus or island)

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
