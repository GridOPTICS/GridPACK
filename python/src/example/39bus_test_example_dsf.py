#!/usr/bin/env python
# -------------------------------------------------------------
# file: 39bus_test_example_dsf.py
# -------------------------------------------------------------
# -------------------------------------------------------------
# Created March 20, 2024 by Perkins
# -------------------------------------------------------------

import sys, os
import gridpack
import numpy as np
from gridpack.dynamic_simulation import DSFullApp, Event, EventVector

# -------------------------------------------------------------
# network_dump_state
# this matches (some) generator power observations output
# -------------------------------------------------------------
def network_dump_state(ds_app):
    nbus = ds_app.totalBuses()
    for bus in range(nbus):
        print(bus,
              ds_app.getBusInfoInt(bus, "BUS_NUMBER"),
              ds_app.getBusInfoString(bus, "BUS_NAME"),
              ds_app.getBusInfoReal(bus, "BUS_VMAG_CURRENT"))
        for g in range(ds_app.numGenerators(bus)):
            print(" gen: ", g,
                  ds_app.getBusInfoString(bus, "GENERATOR_MODEL", g),
                  ds_app.getBusInfoReal(bus, "GENERATOR_PG_CURRENT", g),
                  ds_app.getBusInfoReal(bus, "GENERATOR_QG_CURRENT", g),
                  )
        for l in range(ds_app.numLoads(bus)):
            print("load: ", l,
                  ds_app.getBusInfoInt(bus, "LOAD_NUMBER", l),
                  ds_app.getBusInfoReal(bus, "LOAD_PL_CURRENT", l),
                  ds_app.getBusInfoReal(bus, "LOAD_QL_CURRENT", l)
                  )

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

ds_app.readGenerators(0);
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

ds_app.updateData()
network_dump_state(ds_app)

isteps = 0
while (not ds_app.isDynSimuDone()):
    print("Time = %4.3f" % (float(isteps)*dt))
    print('Before step')
    ds_app.executeOneSimuStep()       
    print('After step')

    if 	isteps%outputob_nsimustep == 0:
        print('Before getObservations')
        ( vMag, vAng,
          rSpd, rAng,
          genP, genQ,
          fOnline, busfreq) = ds_app.getObservations_withBusFreq()
        ob_vals = [ float(isteps)*dt ]
        ob_vals.extend(rSpd)
        ob_vals.extend(rAng)
        ob_vals.extend(genP)
        ob_vals.extend(genQ)
        ob_vals.extend(vMag)
        ob_vals.extend(vAng)
        ob_vals.extend(fOnline)
        ob_vals.extend(busfreq)
        ds_app.updateData()
        network_dump_state(ds_app)
        print('After getObservations')
        
        print('Before insert')
        observation_list.append(ob_vals)
        print('After insert')
    isteps = isteps + 1

network_dump_state(ds_app)
                   
np_data = np.array(observation_list)

import pandas as pd
pd_data = pd.DataFrame(np_data,columns=csvhead)

csv_wrt_f = '39bus_test_example_dsf_observations'

pd_data.to_csv(csv_wrt_f+'.csv', index=False)

ds_app = None
env = None
