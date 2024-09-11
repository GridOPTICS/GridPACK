#!/usr/bin/env python
# -------------------------------------------------------------
# file: dsf.py
# -------------------------------------------------------------
# -------------------------------------------------------------
# Created March 21, 2024 by Perkins
# -------------------------------------------------------------

import sys, os
import gridpack
from gridpack.dynamic_simulation import DSFullApp

# -------------------------------------------------------------
# network_analytics_dump
# -------------------------------------------------------------
def network_analytics_dump(ds_app):
    nbus = ds_app.totalBuses()
    for bus in range(nbus):
        print(bus,
              ds_app.getBusInfoInt(bus, "BUS_NUMBER"),
              ds_app.getBusInfoString(bus, "BUS_NAME"),
              ds_app.getBusInfoInt(bus, "BUS_TYPE"),
              ds_app.numGenerators(bus),
              ds_app.numLoads(bus),
              ds_app.getBusInfoReal(bus, "BUS_VOLTAGE_MAG"))
        for g in range(ds_app.numGenerators(bus)):
            print(" gen: ", g,
                  ds_app.getBusInfoInt(bus, "GENERATOR_NUMBER", g),
                  ds_app.getBusInfoString(bus, "GENERATOR_ID", g),
                  ds_app.getBusInfoReal(bus, "GENERATOR_PG", g),
                  ds_app.getBusInfoReal(bus, "GENERATOR_QG", g),
                  ds_app.getBusInfoReal(bus, "GENERATOR_PG_CURRENT", g),
                  ds_app.getBusInfoReal(bus, "GENERATOR_QG_CURRENT", g),
                  )
        for l in range(ds_app.numLoads(bus)):
            print("load: ", l,
                  ds_app.getBusInfoInt(bus, "LOAD_NUMBER", l),
                  ds_app.getBusInfoString(bus, "LOAD_ID", l),
                  ds_app.getBusInfoReal(bus, "LOAD_PL", l),
                  ds_app.getBusInfoReal(bus, "LOAD_QL", l),
                  ds_app.getBusInfoReal(bus, "LOAD_PL_CURRENT", l),
                  ds_app.getBusInfoReal(bus, "LOAD_QL_CURRENT", l)
                  )
    nbranch = ds_app.totalBranches()
    for branch in range(0, nbranch):
        (f, t) = ds_app.getBranchEndpoints(branch)
        print(branch, f, t, 
              ds_app.getBranchInfoInt(branch, "BRANCH_ELEMENTS"),
              ds_app.getBranchInfoInt(branch, "BRANCH_INDEX"),
              ds_app.getBranchInfoReal(branch, 'BRANCH_FROM_P_CURRENT'),
              ds_app.getBranchInfoReal(branch, 'BRANCH_TO_P_CURRENT'),
              ds_app.getBranchInfoReal(branch, 'BRANCH_FROM_Q_CURRENT'),
              ds_app.getBranchInfoReal(branch, 'BRANCH_TO_Q_CURRENT'))

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

timer = gridpack.CoarseTimer()
t_total = timer.createCategory("Dynamic Simulation: Total Application")
timer.start(t_total)

np = gridpack.NoPrint()
sys.stdout.write("%d: NoPrint status: %r\n" % (comm.rank(), np.status()))
np.setStatus (True)
        
ds_app = DSFullApp()
ds_app.solvePowerFlowBeforeDynSimu(inname, -1)
ds_app.readGenerators();
ds_app.readSequenceData();
ds_app.initialize();
ds_app.setGeneratorWatch();

# Remember the input file was read into the Configuration singleton
conf = gridpack.Configuration()
cursor = conf.getCursor("Configuration.Dynamic_simulation")

faults = ds_app.getEvents(cursor)

ds_app.solvePreInitialize(faults[0])

while (not ds_app.isDynSimuDone()):
    ds_app.executeOneSimuStep()

ds_app.updateData()
network_analytics_dump(ds_app)

timer.stop(t_total)
timer.dump()

# try to force deletion order to avoid problems
del ds_app
del comm
del env

