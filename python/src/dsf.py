#!/usr/bin/env python
# -------------------------------------------------------------
# file: dsf.py
# -------------------------------------------------------------

import sys, os
import gridpack
from gridpack.dynamic_simulation import DSFullApp

# -------------------------------------------------------------
# variable initialization
# -------------------------------------------------------------
program = os.path.basename(sys.argv[0])
usage = "usage: " + program + "input.xml\n"

# -------------------------------------------------------------
# handle command line
# -------------------------------------------------------------
if (len(sys.argv) < 2):
    sys.stderr.write(usage)
    exit(3)

inname = sys.argv[1]

# -------------------------------------------------------------
# main program
# -------------------------------------------------------------

env = gridpack.Environment()

comm = gridpack.Communicator()

np = gridpack.NoPrint()
sys.stdout.write("%d: NoPrint status: %r\n" % (comm.rank(), np.status()))
np.setStatus (True)
        
ds_app = DSFullApp()
ds_app.solvePowerFlowBeforeDynSimu(inname)
ds_app.readGenerators(-1);
ds_app.readSequenceData();
ds_app.initialize();
ds_app.setGeneratorWatch();
ds_app.setup();
ds_app.run();


