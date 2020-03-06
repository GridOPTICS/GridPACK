#!/usr/bin/env python
# -------------------------------------------------------------
# file: hadrec.py
# -------------------------------------------------------------
# -------------------------------------------------------------
# -------------------------------------------------------------
# -------------------------------------------------------------
# Created February 17, 2020 by Perkins
# Last Change: 2020-03-06 06:48:05 d3g096
# -------------------------------------------------------------

import sys, os
import gridpack
import gridpack.hadrec

# -------------------------------------------------------------
# variable initialization
# -------------------------------------------------------------
program = os.path.basename(sys.argv[0])
usage = "usage: " + program + " input.xml\n"

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

hadapp = gridpack.hadrec.Module()
hadapp.solvePowerFlowBeforeDynSimu(sys.argv)
hadapp.transferPFtoDS()
hadapp.initializeDynSimu()

loadshedact = gridpack.hadrec.Action()
loadshedact.actiontype = 0;
loadshedact.bus_number = 5;
loadshedact.componentID = "1";
loadshedact.percentage = -0.2;

(obs_genBus, obs_genIDs, obs_loadBuses, obs_loadIDs, obs_busIDs) = hadapp.getObservationLists()

print (obs_genBus)
print (obs_genIDs)
print (obs_loadBuses)
print (obs_loadIDs)
print (obs_busIDs)

bApplyAct = True
isteps = 0

while (not hadapp.isDynSimuDone()):
    if (bApplyAct and
        ( isteps == 2500 or isteps == 3000 or
          isteps == 3500 or isteps == 4000 )):
        hadapp.applyAction(loadshedact)
    hadapp.executeDynSimuOneStep()
    ob_vals = hadapp.getObservations();
    print (isteps, ob_vals)
    isteps = isteps + 1

btest_2dynasimu = True

if (btest_2dynasimu):
    hadapp.transferPFtoDS()
    hadapp.initializeDynSimu()

    while (not hadapp.isDynSimuDone()):
        if (bApplyAct and
            ( isteps == 2500 or isteps == 3000 or
              isteps == 3500 or isteps == 4000 )):
            hadapp.applyAction(loadshedact)
        hadapp.executeDynSimuOneStep()
        ob_vals = hadapp.getObservations();
        print (isteps, ob_vals)
        isteps = isteps + 1

# It's important to force the deallocation order here
hadapp = None
env = None
