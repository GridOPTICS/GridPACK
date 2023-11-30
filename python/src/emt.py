#!/usr/bin/env python
# -------------------------------------------------------------
# file: emt.py
# -------------------------------------------------------------
# -------------------------------------------------------------
# Created November 29, 2023 by Perkins
# Last Change: 2023-11-30 11:31:08 d3g096
# -------------------------------------------------------------

import sys, os
import gridpack
import gridpack.emt

# -------------------------------------------------------------
# variable initialization
# -------------------------------------------------------------
program = os.path.basename(sys.argv[0])
usage = "usage: " + program + " input.xml"

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

emt = gridpack.emt.EMT()

emt.setconfigurationfile(inname)
emt.solvepowerflow()
emt.setup()

print("start solving: ")
emt.solve()

emt = None
env = None




