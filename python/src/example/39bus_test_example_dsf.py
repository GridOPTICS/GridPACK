#!/usr/bin/env python
# -------------------------------------------------------------
# file: 39bus_test_example_dsf.py
# -------------------------------------------------------------
# -------------------------------------------------------------
# Created March 20, 2024 by Perkins
# -------------------------------------------------------------

import sys, os
import gridpack
from gridpack.dynamic_simulation import DSFullApp, Event, EventVector

# -------------------------------------------------------------
# variable initialization
# -------------------------------------------------------------
program = os.path.basename(sys.argv[0])
usage = "usage: " + program

inname = "input_39bus_step005_v33.xml"

# -------------------------------------------------------------
# main program
# -------------------------------------------------------------
gridpack.NoPrint().setStatus(False)
env = gridpack.Environment()
comm = gridpack.Communicator()

ds_app = DSFullApp()

ds_app.solvePowerFlowBeforeDynSimu(inname, 0)

conf = gridpack.Configuration()
cursor = conf.getCursor("Configuration.Dynamic_simulation")

ds_app.readGenerators();
ds_app.readSequenceData();
ds_app.initialize();
ds_app.setGeneratorWatch();
ds_app.setObservations(cursor)

busfault = Event()

busfault.start = 1.0 # fault start time
busfault.end = 1.1   # fault end time
busfault.step = 0.005  # fault duration simu time step
busfault.isBus = True
busfault.bus_idx = 22 # bus number of the fault

busfaultlist = EventVector([busfault])

(obs_genBus, obs_genIDs, obs_loadBuses, obs_loadIDs, obs_busIDs) = ds_app.getObservationLists()

print (obs_genBus)
print (obs_genIDs)
print (obs_loadBuses)
print (obs_loadIDs)
print (obs_busIDs)

# create observation names for csv file header writting purpose
csvhead = []
csvhead.append('time')
for itmp in range(len(obs_genBus)):
    csvhead.append('gen-spd-%d-%s'%(obs_genBus[itmp], obs_genIDs[itmp]))
for itmp in range(len(obs_genBus)):
    csvhead.append('gen-ang-%d-%s'%(obs_genBus[itmp], obs_genIDs[itmp]))
for itmp in range(len(obs_genBus)):
    csvhead.append('gen-P-%d-%s'%(obs_genBus[itmp], obs_genIDs[itmp]))
for itmp in range(len(obs_genBus)):
    csvhead.append('gen-Q-%d-%s'%(obs_genBus[itmp], obs_genIDs[itmp]))
for itmp in range(len(obs_busIDs)):
    csvhead.append('bus-%d-v'%(obs_busIDs[itmp]))
for itmp in range(len(obs_busIDs)):
    csvhead.append('bus-%d-ang'%(obs_busIDs[itmp]))

print(csvhead)

ds_app.solvePreInitialize(busfaultlist[0]);

observation_list = []

dt = ds_app.getTimeStep()

# output obsevations time step
outputob_time_step = 0.005
outputob_nsimustep = int(outputob_time_step/dt);


isteps = 0
while (not ds_app.isDynSimuDone()):
    print("Time = %4.3f" % (float(isteps)*dt))
    print('Before step')
    ds_app.executeOneSimuStep()       
    print('After step')

    print('Before getObservations')
    ( vMag, vAng,
      rSpd, rAng,
      genP, genQ,
      fOnline, busfreq) = ds_app.getObservations_withBusFreq()
    obs_vals = [ float(isteps)*dt ]
    obs_vals.extend(rSpd)
    obs_vals.extend(rAng)
    obs_vals.extend(genP)
    obs_vals.extend(genQ)
    obs_vals.extend(vMag)
    obs_vals.extend(vAng)
    obs_vals.extend(fOnline)
    obs_vals.extend(busfreq)
    print('After getObservations')
    isteps = isteps + 1

ds_app = None

env = None
