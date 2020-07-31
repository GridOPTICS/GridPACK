#!/usr/bin/env python
# -------------------------------------------------------------
# file: hadrec.py
# -------------------------------------------------------------
# -------------------------------------------------------------
# -------------------------------------------------------------
# -------------------------------------------------------------
# Created February 17, 2020 by Perkins
# Last Change: 2020-07-31 11:37:27 d3g096
# -------------------------------------------------------------

import sys, os
import gridpack
import gridpack.hadrec
import gridpack.dynamic_simulation

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

comm = gridpack.Communicator()

np = gridpack.NoPrint()
sys.stdout.write("%d: NoPrint status: %r\n" % (comm.rank(), np.status()))
# np.setStatus(True)
sys.stdout.write("%d: NoPrint status: %r\n" % (comm.rank(), np.status()))
        
hadapp = gridpack.hadrec.Module()
hadapp.solvePowerFlowBeforeDynSimu(inname, -1)
hadapp.transferPFtoDS()

busfaultlist = gridpack.dynamic_simulation.EventVector()

hadapp.initializeDynSimu(busfaultlist)

if (not np.status()):
    print hadapp.getBusTotalLoadPower(5);
    print hadapp.getBusTotalLoadPower(7);
    print hadapp.getBusTotalLoadPower(9);
    print hadapp.getBusTotalLoadPower(-1); # should be None

    print hadapp.getGeneratorPower(1, "1");
    print hadapp.getGeneratorPower(2, "1");
    print hadapp.getGeneratorPower(3, "1");
    print hadapp.getGeneratorPower(3, "14"); # should be None
    print hadapp.getGeneratorPower(-1, "1"); # should be None
        
loadshedact = gridpack.hadrec.Action()
loadshedact.actiontype = 0;
loadshedact.bus_number = 5;
loadshedact.componentID = "1";
loadshedact.percentage = -0.2;

loadshedact1 = gridpack.hadrec.Action()
loadshedact1.actiontype = 0;
loadshedact1.bus_number = 7;
loadshedact1.componentID = "1";
loadshedact1.percentage = -0.2;

(obs_genBus, obs_genIDs, obs_loadBuses, obs_loadIDs, obs_busIDs) = hadapp.getObservationLists()

if (not np.status()):
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
        hadapp.applyAction(loadshedact1)
    hadapp.executeDynSimuOneStep()
    ob_vals = hadapp.getObservations();
    if (not np.status()):
        print (isteps, ob_vals)
    isteps = isteps + 1

btest_2dynasimu = True

if (btest_2dynasimu):

    busfault = gridpack.dynamic_simulation.Event()
    busfault.start = 10.0;
    busfault.end = 10.2;
    busfault.step = 0.005;
    busfault.isBus = True;
    busfault.bus_idx = 7;

    busfaultlist = gridpack.dynamic_simulation.EventVector([busfault])
    
    hadapp.transferPFtoDS()
    hadapp.initializeDynSimu(busfaultlist)

    while (not hadapp.isDynSimuDone()):
        if (bApplyAct and
            ( isteps == 2500 or isteps == 3000 or
              isteps == 3500 or isteps == 4000 )):
            hadapp.applyAction(loadshedact)
        hadapp.executeDynSimuOneStep()
        ob_vals = hadapp.getObservations();
        if (not np.status()):
            print (isteps, ob_vals)
        isteps = isteps + 1

# See if we could do it again, if we wanted to
hadapp = None
hadapp = gridpack.hadrec.Module()
hadapp.solvePowerFlowBeforeDynSimu(inname)

# It's important to force deallocation order here
hadapp = None
env = None

hadapp = gridpack.hadrec.Module()

 
