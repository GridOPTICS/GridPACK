#!/usr/bin/env python
# -------------------------------------------------------------
# file: pf.py
# -------------------------------------------------------------
# -------------------------------------------------------------
#
# Copyright (c) 2013 Battelle Memorial Institute
# Licensed under modified BSD License. A copy of this license can be found
# in the LICENSE file in the top level directory of this distribution.
#
# -------------------------------------------------------------
# -------------------------------------------------------------
# Created February 10, 2025 by Perkins
# -------------------------------------------------------------

import sys, os
import gridpack
import gridpack.powerflow

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
#
# This mimics pf.x
# -------------------------------------------------------------

env = gridpack.Environment()
comm = gridpack.Communicator()
timer = gridpack.CoarseTimer()

# read and extract some options from the input file
config = gridpack.Configuration()
config.open(inname, comm)
cursor = config.getCursor("Configuration.Powerflow")
useNonLinear = cursor.get("UseNonLinear")
exportPSSE23 = cursor.get("exportPSSE_v23")
exportPSSE33 = cursor.get("exportPSSE_v33")
exportPSSE34 = cursor.get("exportPSSE_v34")
noPrint = cursor.get("suppressOutput")

# Create and initialize power flow application
pfapp = gridpack.powerflow.Powerflow()

if (noPrint):
    pfapp.suppressOutput(True)

pfapp.readNetwork(config, -1)
pfapp.initialize()

if useNonLinear:
    pfapp.nl_solve();
else:
    pfapp.solve();

pfapp.write();
pfapp.saveData;

if (exportPSSE23):
    pfapp.exportPSSE23(exportPSSE23)
if (exportPSSE33):
    pfapp.exportPSSE33(exportPSSE33)
if (exportPSSE34):
    pfapp.exportPSSE34(exportPSSE34)


if not noPrint:
    timer.dump()

# try to force deletion order to avoid problems
del pfapp
del config
del comm
del env
