Real-time path rating fundamental algorithm
=====
The goal of this application is to determine the real-time line rating for the tie-lines
for an area. The real-time line rating is determined through steady-state and dynamic analysis and will yield a maximum MVA flow on the tie-lines without causing any steady-state security
(line flow violations, voltage violations) or dynamic-security (generartor instability) violations. Before, we include 
all sorts of possible variations and violations, let us look at this simple reference version that ONLY considers one tie-line between areas A and B. The only violation considered is the MVA flow violation on tie-line. Here are the steps. (Assume the tie-line is provided by the user)

1. Run base case (assume there is no tie-line flow violation)
2. Run N-1 steady-state contingency analysis removing elements from area A only.
3. Check tie-line flow for each contingency. <b> Note:</b> RATE_B (short-term) line rating which
is typically higher than RATE_A (long-term) rating for the line should be used for
checking violation. The short-term rating allows the line to be loaded higher for shorter period of time.
4. If no violation in step 3, then increase generation in area A and load in area B by 5\%. Go to step 1.
5. Else, there is a tie-line flow violation for some contingency. Mark the generation in area A, load in area B, and the tie line flow for that contingency as the steady-state real-time path rating. If there are multiple contingencies that cause tie-line flows to go over, then mark the one that has the smallest tie-line flow as the steady-state real-time path rating.
6. With the generation, load, and tie-line flow from 5, run the N-1 dynamics simulation.
7. If instability for any contingency, reduce the generation in area A, load in area B by 5% and re-run the
   the N-1 dynamics simulation. The condition for instability is as follows: 
   
   If any generator frequency, 60*(1+ dw), is less than 59 Hz for 0.5 sec and the rate of change of frequency is always negative (i.e., the frequency is dropping), then we consider the system as unstable.
8. If no instability, then we are done. Mark the tie-line flow as the real-time path rating.
Other violations and complex strategies can be incorporated once this basic algorithm is working.


Assume a multi-area system that contains tie-lines between different areas. The
tie-lines are specified by the user and are part of the input data for the
calculation. The application does not determine which lines are tie-lines.

Step 1
======
Run base case power flow calculation. Assume there are no violations. If
violations occur, disable checking on any elements where there are violations.
This will hold for all subsequent calculations

Step 2
======
Run N-1 contingency for the entire system and determine which contingencies cause
a violation.
The list of contingencies will be used to construct the calculations
in step 3.

<b> Note:</b> RATE_B (short-term) line rating which
is typically higher than RATE_A (long-term) rating for the line should be used for
checking violation. The short-term rating allows the line to be loaded higher for shorter period of time.

Step 3
======
If violation exists, scale power generation or load to eleminate the violation for each contingency
that is in violation.
- The contingency represents loss of a line and there is a line violation
  somewhere else in the system. Scale generation downward until violation
  disappears. We may also need to either scale generators and/or loads in the other region up or down to 
  keep power balanced.
- The continency represents loss of a line and there is a voltage violation
  somewhere else in the system. If the voltage is too high, scale generation
  downward until the violation disappears. If the voltage is too low, scale load
  downward until violation disappears.(?)
- The contingency represents loss of a generator and there is a violation of
  some type and we do what??? Same for line contingency.

The generations and load are scaled from between 0 and 100\% of their original
values. Start by an initial scaling of 10\% of all values for each violation to
give an initial set consisting of 10\*Q calculations, where Q is the number of
contingencies that cause a violation. After the initial calculations are run,
run a second series of calculations at increments of 1\% to refine the estimates
from the first set of calculations. This will provide a second set of 10\*Q
calculations. At the end of these calculations, there is a scaling parameter
associated with each contingency that is in violation.

Questions
- Do we take into account area? If the contingency and the violation both occur
  in the same area, should we restrict the scaling to just elements in those
  areas? Should we restrict the scaling to a single area, even if the
  contingency and the violation are in different areas?
  *Yousu*: I think we can pre-define the areas where the adjustable generators
  are located. We can start with the tie-line violations only.
- What do we do if there are multiple violations for a single contingency? How
  do we determine what to scale?
  *Yousu*: The scale will be based on the difference between the generator
  capacity (Pmax or Pmin) and the current generation (Pg). We can scale with a
  weight of Pmax\_i-Pg\_i or uniformly adjust Pg.

Step 4
======
We look at each of the tie-lines and scale the generation at one end and the
load at the other until there is a violation. Then we do what??

*Yousu*: Record the current power flow, especially the changes of
generators/loads. Load shedding is the last step to do.

Step 5
======
We do something with dynamic simulation.

*Yousu*: For the saved case, run DSA to check their transient stability.
