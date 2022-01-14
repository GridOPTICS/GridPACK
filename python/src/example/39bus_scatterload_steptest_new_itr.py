#!/usr/bin/env python
# -------------------------------------------------------------
# file: 39bus_scatterload_steptest_new_itr.py
# -------------------------------------------------------------
# This file shows how to change load during the dynamic simulation
# with GridPACK Hadrec Module
# -------------------------------------------------------------
# -------------------------------------------------------------
# Created July 29, 2021 by Renke Huang
# Last Change: 2020-03-09 12:24:10 d3g096
# -------------------------------------------------------------

import sys, os
import gridpack
import time
import gridpack.hadrec
import numpy as np
import xmltodict
import math
import pandas as pd

pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)

# -------------------------------------------------------------
# variable initialization
# -------------------------------------------------------------

def xml2dict(input_xmlfile):
    file_object = open(input_xmlfile,encoding = 'utf-8')
    try:
        all_the_xmlStr = file_object.read()
    finally:
        file_object.close()
        
    #xml To dict
    convertedDict = xmltodict.parse(all_the_xmlStr)

    return convertedDict
	
inname = "input_39bus_step005_v33_itr.xml"

input_xml_dict = xml2dict(inname)
simu_time_step = float(input_xml_dict["Configuration"]["Dynamic_simulation"]["timeStep"])

# determine the time interval to output the observations
outputob_time_step = 0.01
print ('simu_time_step: ', simu_time_step)
print ('outputob_time_step: ', outputob_time_step)

outputob_nsimustep = int(outputob_time_step/simu_time_step);

# -------------------------------------------------------------
# main program
# -------------------------------------------------------------

print ('before gridpack ini')
env = gridpack.Environment()

hadapp = gridpack.hadrec.Module()

print ('after gridpack ini')

gridpack_starttime = time.time()

hadapp.solvePowerFlowBeforeDynSimu(inname, 0)  # solve the first raw file for power flow, each bus only has one load， 
                                               #the second parameter means the first raw file in xml input file
print ('after gridpack solve pf')
hadapp.transferPFtoDS()

# here we define the fault,
# since this file we just want to show load change,
# we put the fault starting time to be large than the total
# simulation time defined in the input xml file
busfault = gridpack.dynamic_simulation.Event()
busfault.start = 71.0
busfault.end = 71.05
busfault.step = 0.005
busfault.isBus = True
busfault.bus_idx = 401		  

busfaultlist = gridpack.dynamic_simulation.EventVector([busfault])

hadapp.initializeDynSimu(busfaultlist, 0)

# get the observations list, the observations are defined in the input xml file
(obs_genBus, obs_genIDs, obs_loadBuses, obs_loadIDs, obs_busIDs, ob_freqBusIDs) = hadapp.getObservationLists_withBusFreq()

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

bApplyAct = False
isteps = 0

observation_list = []
total_dataconv_time = 0.0

timeser = []

ipd_idx = 1
while (not hadapp.isDynSimuDone()):
        
    hadapp.executeDynSimuOneStep() # execute one time step dynamic simulation
   
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
        # the scatterInjectionLoad will change the loads at the buslist to be the values in the plist and qlist
        #hadapp.scatterInjectionLoadNew(buslist, plist, qlist)
    
   
        bustmp = 508
        ptmp = 800.0
        qtmp = 300.0
        ppu_tmp = ptmp/100.0
        qpu_tmp = qtmp/100.0
        buslist.append(bustmp)
        plist.append(ppu_tmp)
        qlist.append(qpu_tmp)
        # the scatterInjectionLoad will change the loads at the buslist to be the values in the plist and qlist
		# the value of load P and Q should be p.u., based on 100 MVA system base
        hadapp.scatterInjectionLoadNew(buslist, plist, qlist)
    
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
        # the scatterInjectionLoad will change the loads at the buslist to be the values in the plist and qlist
        #hadapp.scatterInjectionLoad(buslist, plist, qlist)
      
        bustmp = 508
        ptmp = 400.0
        qtmp = 100.0
        ppu_tmp = ptmp/100.0
        qpu_tmp = qtmp/100.0
        buslist.append(bustmp)
        plist.append(ppu_tmp)
        qlist.append(qpu_tmp)
        # the scatterInjectionLoad will change the loads at the buslist to be the values in the plist and qlist
        hadapp.scatterInjectionLoadNew(buslist, plist, qlist)
		
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
        # the scatterInjectionLoad will change the loads at the buslist to be the values in the plist and qlist
        #hadapp.scatterInjectionLoad(buslist, plist, qlist)
      
        bustmp = 508
        ptmp = 400.0*1.5
        qtmp = 100.0*1.5
        ppu_tmp = ptmp/100.0
        qpu_tmp = qtmp/100.0
        buslist.append(bustmp)
        plist.append(ppu_tmp)
        qlist.append(qpu_tmp)
        # the scatterInjectionLoad will change the loads at the buslist to be the values in the plist and qlist
        hadapp.scatterInjectionLoadNew(buslist, plist, qlist)
		
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
        # the scatterInjectionLoad will change the loads at the buslist to be the values in the plist and qlist
        #hadapp.scatterInjectionLoad(buslist, plist, qlist)
      
        bustmp = 508
        ptmp = 400.0
        qtmp = 100.0
        ppu_tmp = ptmp/100.0
        qpu_tmp = qtmp/100.0
        buslist.append(bustmp)
        plist.append(ppu_tmp)
        qlist.append(qpu_tmp)
        # the scatterInjectionLoad will change the loads at the buslist to be the values in the plist and qlist
        hadapp.scatterInjectionLoadNew(buslist, plist, qlist)
        
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
        # the scatterInjectionLoad will change the loads at the buslist to be the values in the plist and qlist
        hadapp.scatterInjectionLoadNew(buslist, plist, qlist) # note the buslist, plist, qlist can be lists, which means you can change multiple loads at different buses
    
      
    if 	isteps%outputob_nsimustep == 0:
        before_getob_time = time.time()	
        ob_vals = hadapp.getObservations() # get the observations
        after_getob_time = time.time()	
        total_dataconv_time += (after_getob_time - before_getob_time)

        ob_vals.insert(0,isteps*simu_time_step)	
        observation_list.append(ob_vals)


    #print (isteps, ob_vals)
	
    isteps = isteps + 1

gridpack_endtime = time.time()
print ("-------------!!! total gridpack time :  ", gridpack_endtime - gridpack_starttime)
#print ("-------------!!! total gridpack powerflow and initialization time :  ", gridpack_initime - gridpack_starttime)
#print ("-------------!!! total gridpack data converstion time :  ", total_dataconv_time)

'''
# output the observations
#pd.DataFrame(np.array(observation_list)).to_csv("tamu500_ob_test.csv")
before_csvwr_time = time.time()
#np.savetxt("certs_gridforming_onlygensal_par_adjust_loadchangeact_ob_freq.csv", np.array(observation_list), delimiter=",")
np.savetxt("39bus_gridforming_test_ob.csv", np.array(observation_list), delimiter=",")
total_csvwr_time = time.time() - before_csvwr_time
print ("-------------!!! total csv write time :  ", total_csvwr_time)
'''

before_csvwr_time = time.time()
np_data = np.array(observation_list)

import pandas as pd
pd_data = pd.DataFrame(np_data,columns=csvhead)
#print(pd_data)
#csv_wrt_f = '2gen_test_exdc1_ob_fault'
csv_wrt_f = inname[6:-4] + '_loadchange_new_3_eachstep'

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

total_csvwr_time = time.time() - before_csvwr_time
print("-------------!!! total csv write time :  ", total_csvwr_time)

# It's important to force the deallocation order here

hadapp = None
env = None
