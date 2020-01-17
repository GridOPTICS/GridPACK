Real-time path rating fundamental algorithm
=====
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
