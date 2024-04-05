#!/usr/bin/env python
# -------------------------------------------------------------
# file: 39bus_scatterload_steptest_new_itr.py
# -------------------------------------------------------------
# This file shows how to change load during the dynamic simulation
# with GridPACK Hadrec Module
# -------------------------------------------------------------
# -------------------------------------------------------------
# Created July 29, 2021 by Renke Huang
# Last Change: 2024-04-05 10:12:20 d3g096
# -------------------------------------------------------------

import sys, os
import gridpack
from gridpack.dynamic_simulation import DSFullApp, Event, EventVector
import numpy as np
import math

# -------------------------------------------------------------
# variable initialization
# -------------------------------------------------------------

outputob_time_step = 0.01
inname = "input_39bus_step005_v33_itr.xml"

# -------------------------------------------------------------
# main program
# -------------------------------------------------------------

env = gridpack.Environment()
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
busfault.start = 71.0
busfault.end = 71.05
busfault.step = 0.005
busfault.isBus = True
busfault.bus_idx = 401		  

busfaultlist = EventVector([busfault])

(obs_genBus, obs_genIDs, obs_loadBuses, obs_loadIDs, obs_busIDs, ob_freqBusIDs) = ds_app.getObservationLists_withBusFreq()

print (obs_genBus)
print (obs_genIDs)
print (obs_loadBuses)
print (obs_loadIDs)
print (obs_busIDs)
print (ob_freqBusIDs)

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
for itmp in range(len(ob_freqBusIDs)):
    csvhead.append('bus-%d-freq'%(ob_freqBusIDs[itmp]))

print(csvhead)

ds_app.solvePreInitialize(busfaultlist[0]);

# output obsevations time step
simu_time_step = ds_app.getTimeStep()
outputob_nsimustep = int(outputob_time_step/simu_time_step);

observation_list = []

isteps = 0

while (not ds_app.isDynSimuDone()):

    ds_app.executeOneSimuStep()       
    
    # make step change for the load at 504 and 508
    if isteps == int(5.0/simu_time_step):
        
        bustmp = 504
        ptmp = 700.0
        qtmp = 250.0
        ppu_tmp = ptmp/100.0
        qpu_tmp = qtmp/100.0
        buslist = []
        plist = []
        qlist = []
        buslist.append(bustmp)
        plist.append(ppu_tmp)
        qlist.append(qpu_tmp)
        #ds_app.scatterInjectionLoadNew(buslist, plist, qlist)
    
   
        bustmp = 508
        ptmp = 800.0
        qtmp = 300.0
        ppu_tmp = ptmp/100.0
        qpu_tmp = qtmp/100.0
        buslist.append(bustmp)
        plist.append(ppu_tmp)
        qlist.append(qpu_tmp)
        ds_app.scatterInjectionLoadNew(buslist, plist, qlist)
                
    # make another step change for the load at 504 and 508   
    if isteps == int(15.0/simu_time_step):
        
        bustmp = 504
        ptmp = 300.0
        qtmp = 125.0
        ppu_tmp = ptmp/100.0
        qpu_tmp = qtmp/100.0
        buslist = []
        plist = []
        qlist = []
        buslist.append(bustmp)
        plist.append(ppu_tmp)
        qlist.append(qpu_tmp)
        #ds_app.scatterInjectionLoad(buslist, plist, qlist)
      
        bustmp = 508
        ptmp = 400.0
        qtmp = 100.0
        ppu_tmp = ptmp/100.0
        qpu_tmp = qtmp/100.0
        buslist.append(bustmp)
        plist.append(ppu_tmp)
        qlist.append(qpu_tmp)
        ds_app.scatterInjectionLoadNew(buslist, plist, qlist)
	
    if isteps == int(25.0/simu_time_step):
        
        bustmp = 504
        ptmp = 300.0*1.5
        qtmp = 100.0*1.5
        ppu_tmp = ptmp/100.0
        qpu_tmp = qtmp/100.0
        buslist = []
        plist = []
        qlist = []
        buslist.append(bustmp)
        plist.append(ppu_tmp)
        qlist.append(qpu_tmp)
        #ds_app.scatterInjectionLoad(buslist, plist, qlist)
      
        bustmp = 508
        ptmp = 400.0*1.5
        qtmp = 100.0*1.5
        ppu_tmp = ptmp/100.0
        qpu_tmp = qtmp/100.0
        buslist.append(bustmp)
        plist.append(ppu_tmp)
        qlist.append(qpu_tmp)
        ds_app.scatterInjectionLoadNew(buslist, plist, qlist)
	
    if isteps == int(35.0/simu_time_step):
        
        bustmp = 504
        ptmp = 300.0
        qtmp = 100.0
        ppu_tmp = ptmp/100.0
        qpu_tmp = qtmp/100.0
        buslist = []
        plist = []
        qlist = []
        buslist.append(bustmp)
        plist.append(ppu_tmp)
        qlist.append(qpu_tmp)
        #ds_app.scatterInjectionLoad(buslist, plist, qlist)
      
        bustmp = 508
        ptmp = 400.0
        qtmp = 100.0
        ppu_tmp = ptmp/100.0
        qpu_tmp = qtmp/100.0
        buslist.append(bustmp)
        plist.append(ppu_tmp)
        qlist.append(qpu_tmp)
        ds_app.scatterInjectionLoadNew(buslist, plist, qlist)
        
    # make continuse sin change for the load at 504 and 508
    nstep_change = int(0.005/simu_time_step)	
    if ((isteps >= int(45.0/simu_time_step)) and (isteps%nstep_change == 0)):
        #change load
        buslist = [504, 508]
        
        freq = 1/2.0
        ptmp = 300/100.0* (1.0 + 0.5*math.sin(2*math.pi*freq*isteps*simu_time_step))  # the value of load P and Q should be p.u., based on 100 MVA system base
        qtmp = 100.0/100.0* (1.0 + 0.2*math.sin(2*math.pi*freq*isteps*simu_time_step))
        #print ('%f, %f'%(ptmp, qtmp))
        ptmp1 = 400.0/100.0* (1.0 + 0.5*math.sin(2*math.pi*freq*isteps*simu_time_step))  # the value of load P and Q should be p.u.
        qtmp1 = 100.0/100.0* (1.0 + 0.3*math.sin(2*math.pi*freq*isteps*simu_time_step))
        plist = [ptmp, ptmp1]
        qlist = [qtmp, qtmp1]
        ds_app.scatterInjectionLoadNew(buslist, plist, qlist) 
        
      
    if 	isteps%outputob_nsimustep == 0:
        ( vMag, vAng,
          rSpd, rAng,
          genP, genQ,
          fOnline, busfreq) = ds_app.getObservations_withBusFreq()
        ob_vals = [ float(isteps)*simu_time_step ]
        ob_vals.extend(rSpd)
        ob_vals.extend(rAng)
        ob_vals.extend(genP)
        ob_vals.extend(genQ)
        ob_vals.extend(vMag)
        ob_vals.extend(vAng)
        ob_vals.extend(fOnline)
        ob_vals.extend(busfreq)

        observation_list.append(ob_vals)
	
    isteps = isteps + 1

np_data = np.array(observation_list)
import pandas as pd
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)

pd_data = pd.DataFrame(np_data,columns=csvhead)
csv_wrt_f = inname[6:-4] + '_loadchange_new_3_eachstep_dsf'

pd_data.to_csv(csv_wrt_f+'.csv', index=False)

'''
colVals = pd_data.columns.values.tolist()
dtLen = range(pd_data.shape[0])
fig, ax = plt.subplots()
for colVal in colVals:
    if colVal.find('mag') >= 0:
        ax.plot(dtLen, pd_data[colVal])
#plt.show()
plt.savefig(csv_wrt_f + '.png')
'''

# It's important to force the deallocation order here

ds_app = None
env = None
