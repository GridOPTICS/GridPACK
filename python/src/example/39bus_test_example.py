#!/usr/bin/env python
#
# 39 bus dynamics simulation application
# Can run in serial and parallel
#

import sys, os
import gridpack
import time
import gridpack.hadrec
import numpy as np
import xmltodict
import math
#import pandas as pd
from mpi4py import MPI

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
inname = "input_39bus_step005_v33.xml"

input_xml_dict = xml2dict(inname)

# get the simulation time step from the input file
simu_time_step = float(input_xml_dict["Configuration"]["Dynamic_simulation"]["timeStep"])

# output obsevations time step
outputob_time_step = 0.005
outputob_nsimustep = int(outputob_time_step/simu_time_step);

# -------------------------------------------------------------
# main program
# -------------------------------------------------------------
# Get World communicator
comm = MPI.COMM_WORLD

print ('before gridpack ini')
noprintflag = gridpack.NoPrint()
noprintflag.setStatus(False) # Setting this to True will disable all printing to stdout

# Create GridPACK environment and pass the communicator to it
env = gridpack.Environment()

# Create hadrec module
hadapp = gridpack.hadrec.Module()

gridpack_starttime = time.time()

# solve the power flow
hadapp.solvePowerFlowBeforeDynSimu(inname, 0)  # 0 inidcates that solves the first raw file for power flow, the xml file supports multiple power flow raw files read in

print ('Finished Solving Power Flow')

# transfer data from power flow network to dynamic simulation network
hadapp.transferPFtoDS()

# define a bus fault
busfault = gridpack.dynamic_simulation.Event()
busfault.start = 1.0 # fault start time
busfault.end = 1.1   # fault end time
busfault.step = 0.005  # fault duration simu time step
busfault.isBus = True
busfault.bus_idx = 22 # bus number of the fault		  

busfaultlist = gridpack.dynamic_simulation.EventVector([busfault])

# initialize the dynamic simulation
hadapp.initializeDynSimu(busfaultlist, 0) # 0 inidcates read in the first dyr dynamic parameter file, the xml file supports multiple dyr files read in

# define a line tripping action
linetripact = gridpack.hadrec.Action()
linetripact.actiontype = 1;
linetripact.brch_from_bus_number = 4;
linetripact.brch_to_bus_number = 5;
linetripact.componentID = "1";

# define a genertor tripping action
gentripact = gridpack.hadrec.Action()
gentripact.actiontype = 2;
gentripact.bus_number = 35;
gentripact.componentID = "1";

# define a load shedding action
loadshedact = gridpack.hadrec.Action()
loadshedact.actiontype = 3;
loadshedact.bus_number = 504;
loadshedact.componentID = "1";
loadshedact.percentage = -300;  #  shed 200MW load at bus 501

# get the observation name list 
(obs_genBus, obs_genIDs, obs_loadBuses, obs_loadIDs, obs_busIDs) = hadapp.getObservationLists()

gridpack_initime = time.time()

# print out what kind of observations are defined in the xml file
print (obs_genBus)
print (obs_genIDs)
print (obs_loadBuses)
print (obs_loadIDs)
print (obs_busIDs)

# print network analytics
print("Number of buses:  %d" % (hadapp.totalBuses()))
print("Number of branches: %d" % (hadapp.totalBranches()))
print("Number of generators: %d" % (hadapp.numGenerators()))
print("Number of loads: %d" % (hadapp.numLoads()))
print("Number of lines: %d" % (hadapp.numLines()))
print("Number of storage units: %d" % (hadapp.numStorage()))
print("Branches connected to bus 1: ", hadapp.getConnectedBranches(1))
print("Buses connected to branch 1: ", hadapp.getBranchEndpoints(1))
print("Lines in branch 1: ", hadapp.numLines(1))
print("Generators connected to bus 1: ", hadapp.numGenerators(1))
print("Loads connected to bus 1: ", hadapp.numLoads(1))
sys.exit(1)

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


bApplyAct = True
isteps = 0

observation_list = []
total_dataconv_time = 0.0

timeser = []

while (not hadapp.isDynSimuDone()):
    print('Time = %4.3f' % (isteps*simu_time_step))
    # apply line tripping at t = 15.0 seconds
    #if (bApplyAct and( isteps == 3000 )):
        #hadapp.applyAction(linetripact)

    # apply load shedding at t = 25.0 seconds
    #if (bApplyAct and( isteps == 5000 )):
        #hadapp.applyAction(loadshedact)
		
    # apply generator tripping at t = 35.0 seconds
    #if (bApplyAct and( isteps == 7000 )):
        #hadapp.applyAction(gentripact)

    print('Before step')
    # execute one simulation time step	
    hadapp.executeDynSimuOneStep()       
    print('After step')
    # check if it is the time to get the observations
    if 	isteps%outputob_nsimustep == 0:
        before_getob_time = time.time()
        print('Before getObservations')
        ob_vals = hadapp.getObservations()
        print('After getObservations')
        
        after_getob_time = time.time()	
        total_dataconv_time += (after_getob_time - before_getob_time)
        print('Before insert')
        ob_vals.insert(0,isteps*simu_time_step)	
        observation_list.append(ob_vals)
        print('After insert')
	
    isteps = isteps + 1

gridpack_endtime = time.time()
print ("-------------!!! total gridpack time :  ", gridpack_endtime - gridpack_starttime)
print ("-------------!!! total gridpack powerflow and initialization time :  ", gridpack_initime - gridpack_starttime)
print ("-------------!!! total gridpack data converstion time :  ", total_dataconv_time)

# write the observations into a csv file
before_csvwr_time = time.time()
np_data = np.array(observation_list)

import pandas as pd
pd_data = pd.DataFrame(np_data,columns=csvhead)

csv_wrt_f = '39bus_test_example_observations'

pd_data.to_csv(csv_wrt_f+'.csv', index=False)


total_csvwr_time = time.time() - before_csvwr_time
print("-------------!!! total csv write time :  ", total_csvwr_time)

# It's important to force the deallocation order here

hadapp = None
env = None
