# Real-time path rating

The real-time path rating application is designed to evaluate path rating limits
for a system. The path rating calculation first determines a provisional path
rating based on a series of powerflow contingency analysis calculations. It then
refines this calculation by verifying that the same set of contingencies are
secure based on dynamic simulation. The details of the algorithm are described
in the section at the end of this write-up. Details on running the calulation
are given in the next section.

## Input
The input deck is of the form

```
<?xml version="1.0" encoding="utf-8"?>
<Configuration>
  <RealTimePathRating>
    <printCalcFiles> false </printCalcFiles>
    <sourceArea> 1 </sourceArea>
    <sourceZone> 1 </sourceZone>
    <destinationArea> 1 </destinationArea>
    <destinationZone> 2 </destinationZone>
    <calculateGeneratorContingencies>true</calculateGeneratorContingencies>
    <calculateLineContingencies>true</calculateLineContingencies>
    <contingencyList>contingencies_14.xml</contingencyList>
    <useBranchRatingB>false</useBranchRatingB>
    <maxVoltage>1.1</maxVoltage>
    <minVoltage>0.9</minVoltage>
    <checkQLimit>false</checkQLimit>
    <monitorGenerators> true </monitorGenerators>
    <frequencyMaximum> 61.8 </frequencyMaximum>
    <contingencyDSStart> 1.0 </contingencyDSStart>
    <contingencyDSEnd> 1.015</contingencyDSEnd>
    <contingencyDSTimeStep> 0.005 </contingencyDSTimeStep>
    <!--
    <tieLines>
      <tieLine>
        <Branch> 2 5 </Branch>
        <Tag> BL </Tag>
      </tieLine>
      <tieLine>
        <Branch> 9 14 </Branch>
        <Tag> BL </Tag>
      </tieLine>
    </tieLines>
    -->
  </RealTimePathRating>
  <Powerflow>
    <networkConfiguration> IEEE14_ca_mod_rate6.raw </networkConfiguration>
    <maxIteration>50</maxIteration>
    <tolerance>1.0e-3</tolerance>
    <LinearSolver>
      <PETScOptions>
        -ksp_type richardson
        -pc_type lu
        -pc_factor_mat_solver_package superlu_dist
        -ksp_max_it 1
      </PETScOptions>
    </LinearSolver>
  </Powerflow>
  <Dynamic_simulation>
    <generatorParameters>IEEE14.dyr</generatorParameters>
    <simulationTime>10.0</simulationTime>
    <timeStep>0.005</timeStep>
    <faultEvents>
      <faultEvent>
        <beginFault> 1.0</beginFault>
        <endFault>   1.03</endFault>
        <faultBranch>6 7</faultBranch>
        <timeStep>   0.005</timeStep>
      </faultEvent>
    </faultEvents>
    <generatorWatch>
      <generator>
       <busID> 1 </busID>
       <generatorID> 1 </generatorID>
      </generator>
    </generatorWatch>
    <generatorWatchFrequency> 1 </generatorWatchFrequency>
    <generatorWatchFileName> gen_watch.csv </generatorWatchFileName>
    <LinearMatrixSolver>
      <!--
                   These options are used if SuperLU was built into PETSc 
      -->
      <Ordering>nd</Ordering>
      <Package>superlu_dist</Package>
      <Iterations>1</Iterations>
      <Fill>5</Fill>
    </LinearMatrixSolver>
  </Dynamic_simulation>
</Configuration>
```

This file contains three major blocks, the `RealTimePathRating` block, the
`Powerflow` block and the `Dynamic_simulation` block. The `Powerflow` and
`Dynamic_simulation` blocks are described elsewhere and are only mentioned here
to note that the PSS\E raw file describing the system is specified in the
`Powerflow` `networkConfiguration` field and the PSS\E dyr file describing the
generator properties is specified in the `Dynamic_simulation`
`generatorParameters` field. The fields in the `RealTimePathRating` block are
described below.

- `printCalcFiles`: This flag indicates whether a file is printed describing
each power flow calculation for each contingency for each value of the
rating parameter. For a reasonable sized system, this will likely result
and enormous number of files so it should be set to false in most cases.

- `sourceArea`,`sourceZone`: Specifies the zone within an area that
represents the generation source in the calculation. If no zone is
specified, then the entire area is assumed to be the source.

- `destinationArea`,`destinationZone`: Specifies the zone within an area that
represents the destination load in the calculation. If no zone is
specified, then the entire area is assumed to be the destination.

- `calculateGeneratorContingencies`, `calculateLineContingencies`:
if true, automatically generate a list of generator and line contingencies
by assuming that all active generators and/or lines contribute a contingency.
If both of these are false or not specified, assume that contingencies are
specified in and external contingency file.

- `contingencyList`: specifies an external contingency file that contains a
listing of all contingencies that are used in the calculation. This file
uses the same format as the contingency application. If either the
`calculateGeneratorContingencies` or `calculateLineContingencies` are set
to true, then this file is ignored.

- `useBranchRatingB`: use the RATEB parameter to evaluate line contingencies
instead of RATEA.

- `minVoltage`,`maxVoltage`: minimum and maximum values of the voltage used
in evaluating power flow contingencies.

- `checkQLimit`: perform the Q-limit test, convert PV
buses to PQ buses if the test fails and rerun power flow calculation.

- `monitorGenerators`: monitor generators for frequency violations. If this
parameter is false then the dynamic simulation part of the path rating
calculation is skipped. Only generators in the source zone or area are
checked.

- `frequencyMaximum`: maximum allowable frequency for monitored generators.
The default is 61.8 Hz.

- `contingencyDSStart`,`contingencyDSEnd`,`contingencyDSTimeStep`: Start
time, end time and time step for faults in the dynamic simulation
contingecies. These apply to all contingencies evaluated in the dynamic
simulation phase of the path rating calculation.

- `tieLines`: This field is used to define user-specified tie-lines. In the
  event that this field is not specified, the RTPR will calculate tie-lines
  automatically by choosing all active lines between the source and destination
  areas. Each tie-line is specified by the `tieLine` sub-block, which contains two
  additional fields. The `Branch` field contains two integers, representing the
  indices of the bus at each end of the line and the `Tag` field contains a one or
  two character string representing the line ID.

## Output

Output for the real-time path rating calculation is fairly compact. If
`printCalcFiles` parameter is set to false, then the main result is the path
rating values in standard output. If standard out is redirected to a file
(recommended) then the rating values can be found by searching for `Rating` in
the output. This will typically yield something like

```
% grep Rating file.out
Final Power Flow Rating: 1.850000
Final Dynamic Simulation Rating: 1.640000
```

If the system hits an upper or lower bound when performing the rating
calculations, some additional lines may appear in the output (these can also be
found by searching on `Rating`).

In addition to the rating values, the calculation exports a set of files with
the names `line_flt_cnt_XXX.txt`, where `XXX` stands for a value of the rating
parameter. These are summaries of each of the power flow contingency
calculations describing how many times each line exhibited a fault for each
value of the rating parameter that was tested in the simulation. Note that the
rating calculation itself only checks tie-lines while these files list faults on
all lines in the system.

## Fundamental algorithm

The goal of this application is to determine the real-time line rating (RTPR) for the
tie-lines between two areas or zones. The real-time line rating is determined through
steady-state and dynamic analysis and will yield a maximum MVA flow on the tie-lines
without causing any steady-state security (line flow violations, voltage violations)
or dynamic-security (generator instability) violations. This implementation considers a
simple model that only considers tie lines between two areas or zones. The only
violation considered is MVA flow violations on tie-lines for the power flow
calculation and frequency violations for the dynamic simulation. The implementation
consists of the following steps:

1. The system is divided into a source region and a destination region. The
destination region is characterized by the loads in that region and the source region
is characterized by the generators in that region. The destination and source regions
can be characterized by either an area or a zone within an area.

2. Tie-lines between the source and destination regions can be specified in one
of two ways. The first is that the user can manually specify each of the lines
that they want to consider as tie-lines. The second method is to have the RTPR
application itself identify the tie-lines. In this case, any active line
directly connecting the source and destination areas will be considered a
tie-line.

3. The RTPR calculation proceeds by determing a scaling factor that is used to
simultaneously increase or decrease the loads in the destination region while
simultaneously matching the change in load with a change in generation in the
source region. The path rating is the maximum value of the scaling parameter
that allows operation with no overloads on the tie-lines and no large deviations
in frequency from the generators.

4. For a given value of the rating parameter, the active loads in the
destination region are scaled by the rating parameter. This will yield an
absolute increase or decrease in the power consumption in the destination
region. The change in power consumption is matched by a corresponding change in
the power generation in the source region. The change in power generation is
evaluated by the following algorithm:
   1. Calculate the total real power generation in the source region for the active
generators for the unperturbed system.
   2. Calculate that margin for the generators in the source region. If the rating
parameter is greater than 1, then calculate the total real power generation
based on the maximum allowable generation. If the rating parameter is less than
1, then calculate the minimum power generation based in the minimum allowable
generation for each generator. The difference between this number and the original power
generation is the margin.
   3. If the change in load is greater than the available margin then stop the
path rating calculation and calculate the value of the rating parameter based on
the available margin. If the change in load is less than the margin, then
calculate the fraction of the margin that will match the change in load and then
change the generation for the individual generators by adding (or subtracting)
the previously determined fraction times the individual margin for the generator
from the unperturbed generation value.

5. The path rating is determined by evaluating the system against a series of
N-1 contingencies. These contingencies can be specified either explicitly by the
user or by having the RTPR application determine them automatically. For
automatic generation, it is possible to choose only generator faults, only line
faults or both. For automatically created generator faults, a list of all active
generators in the source regions is determined and each of these generators is
associated with a fault. For line faults, all active lines connecting two buses
that are both located in the source area or zone are associated with a fault.

6. A complete set of power flow calculations containing both the base case (no
contingency) and the complete set of contingencies is run at a given value of
the rating parameter, starting with a rating of 1.0. For both the base case and
each contingency, the tie lines are checked to see if there is an overflow
violation. The overflow violation can be based on the Rating A or Rating B
parameter. The default for the RTPR calculation is Rating B. If there is a
violation for any contingency or the base case, then the rating parameter is
decreased until there are no violations. If there is no violation initially, then
the rating parameter is increased until there is a violation on the tie-lines.
The final value of the rating parameter is the power flow rating for the system.
The calculation proceeds by changing the path rating in units of 0.05 (5\%). When
a limit is found, the calculation reverts to the previous value of the rating and
start changing in increments of 0.01 (1\%). The use of two different increments
is designed to speed up the calculation for final values of the rating parameter
that are significantly different from the starting value of 1.0.

7. After the power flow rating is completed, the system is then checked using
dynamical simulation. The dynamic simulation rating is assumed to be more
stringent than the power flow rating, so the power flow rating is considered to
be an upper bound on the dynamic simulation value. The dynamic simulation rating
will only check for lower values of the rating.

8. Each of the contingencies created for the power flow rating is converted to a
fault for dynamic simulation. These faults must be supplemented by information
on when the fault occurs and how long it lasts. The dynamic simulation path
rating defines a start time for the fault, an end time for the fault and a
simulation length that is applied to all faults. Simulations are run for each of
the faults and the frequency of all generators in the system is checked to see
if any frequencies exceed a threshold. This threshold is configurable but is
currently set to 61.8 Hz. If there is a frequency violation for any fault, the
rating parameter is considered to generate a violation.

9. The first dynamic simulation path rating calculation starts with the power
flow rating value. If there are no frequency violations, then the calculation
stops and the dynamic simulation path rating is the same as the power flow path
rating. If there is a violation, then the rating is decreased until there is no
violation. This value is then the dynamic simulation path rating. The rating
parameter for dynamic simulation is used to scale the source and destination
areas or zones in the same way as for the power flow ratings.
