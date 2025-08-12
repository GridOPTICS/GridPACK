#!/usr/bin/env python
# -------------------------------------------------------------
# file: stes.py
# -------------------------------------------------------------
# -------------------------------------------------------------
# Created August 11, 2025 by Perkins
# Last Change: 2025-08-12 08:21:44 d3g096
# -------------------------------------------------------------

import sys, os
from optparse import OptionParser
import gridpack
import gridpack.state_estimation

# -------------------------------------------------------------
# variable initialization
# -------------------------------------------------------------
program = os.path.basename(sys.argv[0])
usage = "usage: " + program

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

timer = gridpack.CoarseTimer()

config = gridpack.Configuration()
config.open(inname, comm)

se = gridpack.state_estimation.SEApp()

se.readNetwork(config)
se.initialize()
se.readMeasurements()
se.solve()
se.saveData()
se.write()

if not se.hasConverged():
    sys.stderr("%s: warning: State estimation did not fully converge, but results were written anyway.\n")

timer.dump()

del se
del config
del comm
del env

