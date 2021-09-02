#!/usr/bin/env python
# -------------------------------------------------------------
# file: 39bus_test_example.py
# -------------------------------------------------------------
# -------------------------------------------------------------
# -------------------------------------------------------------
# -------------------------------------------------------------
# Created June 4, 2021 by Renke Huang
# -------------------------------------------------------------
# -------------------------------------------------------------

import sys, os
import gridpack
import time
import gridpack.hadrec
import numpy as np
import xmltodict
import math
#import pandas as pd

def xml2dict(input_xmlfile):
    file_object = open(input_xmlfile,encoding = 'utf-8')
    try:
        all_the_xmlStr = file_object.read()
    finally:
        file_object.close()
        
    #xml To dict
    convertedDict = xmltodict.parse(all_the_xmlStr)

    return convertedDict
	
# define input file for gridpack hadrec simulator
inname = "input_39bus_step005_pss.xml"

input_xml_dict = xml2dict(inname)

# get the simulation time step from the input file
simu_time_step = float(input_xml_dict["Configuration"]["Dynamic_simulation"]["timeStep"])

# set up the pss bus number and genid
pssbusno = 34
genid = '1'

#set up wide area bus numbers for get frequency
freqbuslist = [34, 30]

# whether use the wide area control signals
useWideAreaControl = True

# output obsevations time step
outputob_time_step = 0.005
outputob_nsimustep = int(outputob_time_step/simu_time_step);

# -------------------------------------------------------------
# main program
# -------------------------------------------------------------

print ('before gridpack ini')
noprintflag = gridpack.NoPrint()
noprintflag.setStatus(True) # if you set this to be True, no any print out from GridPACK Hadrec module, usefull for parallel running with ray to disable any temp print outs
run_gridpack_necessary_env = gridpack.Environment()
hadapp = gridpack.hadrec.Module()

print ('after gridpack ini')

gridpack_starttime = time.time()

# solve the power flow
hadapp.solvePowerFlowBeforeDynSimu(inname, 0)  # 0 inidcates that solves the first raw file for power flow, the xml file supports multiple power flow raw files read in

print ('after gridpack solve pf')

# transfer data from power flow network to dynamic simulation network
hadapp.transferPFtoDS()

# define a bus fault
busfault = gridpack.dynamic_simulation.Event()
busfault.start = 1.0 # fault start time
busfault.end = 1.1   # fault end time
busfault.step = 0.005  # fault duration simu time step
busfault.isBus = True
busfault.bus_idx = 16 # bus number of the fault		  

busfaultlist = gridpack.dynamic_simulation.EventVector([busfault])

# initialize the dynamic simulation
hadapp.initializeDynSimu(busfaultlist, 0) # 0 inidcates read in the first dyr dynamic parameter file, the xml file supports multiple dyr files read in

# get the observation name list 
(obs_genBus, obs_genIDs, obs_loadBuses, obs_loadIDs, obs_busIDs, ob_freqBusIDs) = hadapp.getObservationLists_withBusFreq()

gridpack_initime = time.time()

# print out what kind of observations are defined in the xml file
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
        
# find the index in the ob array for the bus frequency
freqbuslist_obidx = []
for busnotmp in freqbuslist:
    freqbusobidxtmp = ob_freqBusIDs.index(busnotmp)
    freqbuslist_obidx.append(freqbusobidxtmp)

print ('freqbuslist_obidx:', freqbuslist_obidx)
#bApplyAct = False
isteps = 0

observation_list = []
total_dataconv_time = 0.0

timeser = []

'''
 the order of the observations are:
    gen spd, size is len(obs_genBus)
    gen ang, size is len(obs_genBus)
    gen P, size is len(obs_genBus)
    gen Q, size is len(obs_genBus)
    bus voltage magnitude, size is len(obs_busIDs)
    bus voltage phase angle, size is len(obs_busIDs)
    bus frequency, size is len(ob_freqBusIDs)
    
'''

ob_vals = hadapp.getObservations()
nobfreqlen = len(ob_freqBusIDs)

while (not hadapp.isDynSimuDone()):
    
    # here we get the bus frequency from bus 34 and bus 30
    wideareafreq_sendtohelics = ob_vals[-(nobfreqlen-freqbuslist_obidx[0])] - ob_vals[-(nobfreqlen-freqbuslist_obidx[-1])]
    
    # Dexing, please add helics code to send out the bus frequency signals to NS3
    
    # Dexing, please add helics code to receive the bus frequency signals from NS3
    
    wideareafreq_recfromhelics = wideareafreq_sendtohelics
    #set the wide area frequency for the pss
    if useWideAreaControl:
        hadapp.setWideAreaControlSignal(pssbusno, genid, wideareafreq_recfromhelics)
        
    #print ('-------python set wide area, wide area: %.6f, freq1: %.6f, freq2: %.6f'%(ob_vals[-2], ob_vals[-1], wideareafreq))
    
    # execute one simulation time step	
    hadapp.executeDynSimuOneStep()       

    # get the observations
    ob_vals = hadapp.getObservations()

    ob_vals.insert(0,isteps*simu_time_step)	
    observation_list.append(ob_vals)

    isteps = isteps + 1

gridpack_endtime = time.time()
print ("-------------!!! total gridpack time :  ", gridpack_endtime - gridpack_starttime)

# write the observations into a csv file
before_csvwr_time = time.time()
np_data = np.array(observation_list)

import pandas as pd
pd_data = pd.DataFrame(np_data,columns=csvhead)

if useWideAreaControl:
    csv_wrt_f = '39bus_test_example_with_wideareacontrolsignal_observations'
else:
    csv_wrt_f = '39bus_test_example_without_wideareacontrolsignal_observations'

pd_data.to_csv(csv_wrt_f+'.csv', index=False)

total_csvwr_time = time.time() - before_csvwr_time
print("-------------!!! total csv write time :  ", total_csvwr_time)

# It's important to force the deallocation order here

hadapp = None
run_gridpack_necessary_env = None
