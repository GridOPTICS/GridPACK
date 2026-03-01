#!/usr/bin/env python
# -------------------------------------------------------------
# file: stes.py
# -------------------------------------------------------------
# -------------------------------------------------------------
# Created August 11, 2025 by Perkins
# Last Change: 2025-09-23 12:06:59 d3g096
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

measname = config.get("Configuration.State_estimation.measurementList")
print(measname)

mconfig = gridpack.Configuration()
mconfig.open(measname, comm)


measures = se.getMeasurements(mconfig)

# List some measurements
for m in measures:
    if m.type() == "VM":
        print("%s: Bus %d: %f, %f" % (m.type(), m.busid, m.value, m.deviation))

se.readNetwork(config)
se.initialize()

#se.readMeasurements()
se.setMeasurements(measures)

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

