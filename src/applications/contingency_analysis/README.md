The contingency analysis calculation produces a number of files summarize the
results of the entire set of individual contingency simulations. Only results
for contingencies that ran to completion are included. Calculations that failed
either because of a numerical instability or because the calculations failed to
converge are not included in the results. The output files are described below.

**success.txt**: This file summarizes that results of each contingency and
reports 1) whether the contingency calculation successfully ran to completion
and 2) whether a violation was found. If a violation is found, the calculation
reports on whether it was a on a bus, on a branch, or both.

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

**pq\_changed\_cnt.txt** This file is only created if the checkQLimit flag is
set to "true" in the input file. It counts the number of times a PV bused is
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
