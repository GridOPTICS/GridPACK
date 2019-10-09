Assume a multi-area system that contains tie-lines between different areas. The
tie-lines are specified by the user and are part of the input data for the
calculation. The application does not determine which lines are tie-lines.

#Step 1
Run base case power flow calculation. Assume there are no violations. If
violations occur, disable checking on any elements where there are violations.
This will hold for all subsequent calculations

#Step 2
Run N-1 contingency for the entire system and determine which contingencies cause
a failure. The list of contingencies will be used to construct the calculations
in step 3. Questions:
- Are we paying any attention to area. For example, d

#Step 3
Scale power generation or load to eleminate the violation for each contingency
that is in violation.
- The contingency represents loss of a line and there is a line violation
  somewhere else in the system. Scale generation downward until violation
  disappears.
- The continency represents loss of a line and there is a voltage violation
  somewhere else in the system. If the voltage is too high, scale generation
  downward until the violation disappears. If the voltage is too low, scale load
  downward until violation disappears.(?)
- The contingency represents loss of a generator and there is ???

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
- What do we do if there are multiple violations for a single contingency? How
  do we determine what to scale?

#Step 4
We look at each of the tie-lines and scale the generation at one end and the
load at the other until there is a violation. Then we do what??

#Step 5
We do something with dynamic simulation.
