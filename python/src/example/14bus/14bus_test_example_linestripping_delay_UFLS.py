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
import cmath
import pandas as pd

# define a function that could read in the input xml contents as dictionary
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
inputfile = "input_14bus_step005_v33.xml"
tielinefile = "14bus_branchdata_tieline.csv"


'''
!!!!!!!!!!!!!!
For the case -1, set ipfcase = 0, and bloadshedding = False
For the case -2, set ipfcase = 4, and bloadshedding = False
For the case -3, set ipfcase = 4, and bloadshedding = True
!!!!!!!!!!!!!!
'''

'''
For the ipfcase variable, we define different power flow raw files 
with different line ratings, for tripping different lines experiments

choose 0, no line trip
choose 2, only trip 2 tie-lines, 5-6 and 4-9
choose 3, trip 3 tie-lines, 5-6, 4-7 and 4-9
choose 4, trip 4 tie-lines, 5-6, 4-9, 10-11, and 13-14
'''

ipfcase = 4

'''
# this variable defines the dyr file index in the xml input file, 
you may add more dyr files for your own testing
'''
idyrcase =0  

'''
# whether implement manual load shedding to mimic the energy storage system 
# increase output power
'''
bloadshedding = True #False


# under frequency load shedding parameters
bUFLS = True
loadshedingbuslist = [2,3,4,5,6,9,10,11,12,13,14]  # the load buses
freqthre = 59.5  # under frequency load shedding threshold
loadshedperc = 0.1  # load shedding percentage
loadshedingtimer = [-1 for i in range (len(loadshedingbuslist))] # timer to count for the load shedding delay
loadshedingtimersteps = 40 # load shedding delay time steps

# line tripping relay setting
blinetriprelay = True #False
linetriptimersteps = 35

input_xml_dict = xml2dict(inputfile)

# get the simulation time step from the input file
simu_time_step = float(input_xml_dict["Configuration"]["Dynamic_simulation"]["timeStep"])

# output observations time step
outputob_time_step = 0.005
outputob_nsimustep = int(outputob_time_step/simu_time_step);

# -------------------------------------------------------------
# main program
# -------------------------------------------------------------

print ('before gridpack ini')
noprintflag = gridpack.NoPrint()
noprintflag.setStatus(False) # if you set this to be True, no any print out from GridPACK Hadrec module
run_gridpack_necessary_env = gridpack.Environment()
hadapp = gridpack.hadrec.Module()

print ('after gridpack ini')

gridpack_starttime = time.time()

'''
# solve the power flow
# ipfcase inidcates the index for choosing which raw file for power flow, 
the xml input file supports multiple power flow raw files read in

'''

hadapp.solvePowerFlowBeforeDynSimu(inputfile, ipfcase)  
print ('after gridpack solve pf')

#----------read in tie-line data information csv input file here--------------
linedf = pd.read_csv(tielinefile)

'''

In the following for loop, do the following three tasks before simulation:
1. Construct the line flow output csv file head
2. read the line parameters from the input data file
3. Compute the Y matrix of each tie-line

The variable lineparalist will store Tieline parameters need for computing the line flow and 
checking whether the line flow exceeds the line MVA rating
This variable list is two dimensional
First dimension: Tie-line index, it is the same tie-line order from the tie-line input file (tielinefile)
Second dimension: parameters for each line:
0: from bus, 1: to bus, 2: ckt id, 3:R,  4:X, 5: B, 6: TAP, 7: Phase Shift, 8: Rating

The Y11, Y12, Y21, Y22 are the four elements for the Y matrix of the single line
[ Y11 Y12
  Y21 Y22]
  
Definition and how to compute the Y matrix of the line from line parameters,
Please refer to the Matpower 5.0 User's Manual, Section 3.2 Branch Model
https://matpower.org/docs/MATPOWER-manual-5.0.pdf

'''
lineparalist = []
Y11 = []
Y12 = []
Y21 = []
Y22 = []

lineflow_csvhead = []
lineflow_csvhead.append('time')

nlinestatus = [1 for i in range(len(linedf))]
linetriptimer = [-1 for i in range (len(linedf))]
for iline in range(len(linedf)):
    singlelinepara = []
    brfrombus  =  linedf['From Bus  Number'][iline]# from bus number
    brtobus = linedf['To Bus  Number'][iline] # to bus number
    cktid = str(linedf['Id'][iline])  # branck ckt ID
    
    lineflow_csvhead.append('branch-flow-P-%d-%d-%s'%(brfrombus, brtobus, cktid))
    lineflow_csvhead.append('branch-flow-Q-%d-%d-%s'%(brfrombus, brtobus, cktid))
    
    singlelinepara.append(brfrombus)
    singlelinepara.append(brtobus)
    singlelinepara.append(cktid)
    
    branchpar = 'BRANCH_R'  # the parameter of the branch we want to get value from the simulator
    testpar = hadapp.getDataCollectionBranchParam(brfrombus, brtobus, cktid, branchpar)
    singlelinepara.append(testpar)
    
    branchpar = 'BRANCH_X'
    testpar = hadapp.getDataCollectionBranchParam(brfrombus, brtobus, cktid, branchpar)
    singlelinepara.append(testpar)
    
    branchpar = 'BRANCH_B'
    testpar = hadapp.getDataCollectionBranchParam(brfrombus, brtobus, cktid, branchpar)
    singlelinepara.append(testpar)
    
    branchpar = 'BRANCH_TAP'
    testpar = hadapp.getDataCollectionBranchParam(brfrombus, brtobus, cktid, branchpar)
    #print ('original branch tap: ', testpar)
    if testpar == 0.0:
        testpar = 1.0
    singlelinepara.append(testpar)
    
    branchpar = 'BRANCH_SHIFT'
    testpar = hadapp.getDataCollectionBranchParam(brfrombus, brtobus, cktid, branchpar)
    singlelinepara.append(testpar)
    
    branchpar = 'BRANCH_RATING_A'
    testpar = hadapp.getDataCollectionBranchParam(brfrombus, brtobus, cktid, branchpar)
    singlelinepara.append(testpar)

    
    '''
    print ('       index %4d, Branch from bus %4d to bus %4d, CKT: %2s, R: %7.5f, X: %7.5f, B: %7.5f, TAP: %5.4f, SHIFT: %7.4f, RATING-_A: %7.2f '%
           (iline, brfrombus, brtobus, cktid, singlelinepara[3],  singlelinepara[4], singlelinepara[5], singlelinepara[6], singlelinepara[7], singlelinepara[8]))
    '''
    
    R = singlelinepara[3]
    X = singlelinepara[4]
    B = singlelinepara[5]
    TAP = singlelinepara[6]
    SHIFT = singlelinepara[7]
    
    lineparalist.append(singlelinepara)
    
    '''
    # Compute  the four elements of the Y matrix from line parameters,
    Please refer to the Matpower 5.0 User's Manual, Section 3.2 Branch Model
    https://matpower.org/docs/MATPOWER-manual-5.0.pdf
    '''

    ys = complex(1.0, 0.0)/complex(R, X)
    jay = complex(0.0, 1.0)
    y11= (ys+jay*B/2.0)/(TAP*TAP)
    y12 = -ys/(TAP*cmath.exp(-jay*SHIFT))
    y21 = -ys/(TAP*cmath.exp(jay*SHIFT))
    y22 = (ys+jay*B/2.0)
    Y11.append(y11)
    Y12.append(y12)
    Y21.append(y21)
    Y22.append(y22)

# get the load original normal operation values of the system
# set the load value output csv file head
loadvalue_csvhead = []
loadvalue_csvhead.append('time')
loadorgvaluelist = []
for ibus in (loadshedingbuslist):  
    loadpar = 'LOAD_PL'
    testpar = hadapp.getDataCollectionLoadParam(ibus, '1', loadpar)
    loadorgvaluelist.append(testpar)
    loadvalue_csvhead.append('bus-%d-load-P'%(ibus))
    
loadremainingvaluelist = loadorgvaluelist[:] # list record the remaining load of each load bus
    
# transfer data from power flow network to dynamic simulation network
hadapp.transferPFtoDS()

'''
Here we show how to define a three-phase bus fault and 
How to tell the simulator when and where the fault should be implemented

Note: in this simulation case, we do not want to apply the bus fault,
so we set the bus fault start time to be longer than the simulation time,
to avoid the bus fault to be applied in the simulation
'''
busfault = gridpack.dynamic_simulation.Event()
busfault.start = 100.0 # fault start time
busfault.end = 100.1   # fault end time
busfault.step = 0.005  # fault duration simu time step
busfault.isBus = True
busfault.bus_idx = 22 # bus number of the fault		  

busfaultlist = gridpack.dynamic_simulation.EventVector([busfault])

'''
# initialize the dynamic simulation
# idyrcase indcates the index for choosing which dyr dynamic file for power flow, 
the xml input file supports multiple dyr dynamic files read in
'''
hadapp.initializeDynSimu(busfaultlist, idyrcase) 


# define a generator trip action
gentripact = gridpack.hadrec.Action()
gentripact.actiontype = 2;
gentripact.bus_number = 8;
gentripact.componentID = "1";

# define a load shedding action
loadshedact = gridpack.hadrec.Action()
loadshedact.actiontype = 3;
loadshedact.bus_number = 11;
loadshedact.componentID = "1";
loadshedact.percentage = -65.0/3.0/4;  

loadshedact2 = gridpack.hadrec.Action()
loadshedact2.actiontype = 3;
loadshedact2.bus_number = 12;
loadshedact2.componentID = "1";
loadshedact2.percentage = -65.0/3.0/4; 

loadshedact3 = gridpack.hadrec.Action()
loadshedact3.actiontype = 3;
loadshedact3.bus_number = 13;
loadshedact3.componentID = "1";
loadshedact3.percentage = -65.0/3.0/4; 

loadshedact4 = gridpack.hadrec.Action()
loadshedact4.actiontype = 3;
loadshedact4.bus_number = 9;
loadshedact4.componentID = "1";
loadshedact4.percentage = -10.0/4; 

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

# form some dictionary to record the index of the bus voltage magnitude and 
# bus voltage phase angle observations in the whole observation array

busmagidx_dict = {}
busangidx_dict = {}
for itmp in range (len(obs_busIDs)):
    busmagidx_dict[obs_busIDs[itmp]] = 4*len(obs_genBus)+itmp
    busangidx_dict[obs_busIDs[itmp]] = 4*len(obs_genBus)+len(obs_busIDs)+itmp
    
# form some dictionary to record the index of the bus frequency
# observations in the whole observation array

busfreqidx_dict = {}
for itmp in range (len(loadshedingbuslist)):
    
    busnumtmp = loadshedingbuslist[itmp] 
    busfreqidx_dict[busnumtmp] = 4*len(obs_genBus)+2*len(obs_busIDs)+ob_freqBusIDs.index(busnumtmp)
    
#print ('busmagidx_dict: ', busmagidx_dict)

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

bApplyGenTripAct = True
isteps = 0

observation_list = []
lineflow_list = []
loadvalue_list = []
total_dataconv_time = 0.0

timeser = []

'''
The following While loop is the main dynamic simulation loop

'''
while (not hadapp.isDynSimuDone()):

    # apply generator tripping at t = 5.0 seconds
    if (bApplyGenTripAct and( isteps == 200 )):
        hadapp.applyAction(gentripact)
        
    # apply load shedding at t = 25.0 seconds
    if (bloadshedding and( isteps == 205 )):
        hadapp.applyAction(loadshedact)
        hadapp.applyAction(loadshedact2)
        hadapp.applyAction(loadshedact3)
        hadapp.applyAction(loadshedact4)
      
    if (bloadshedding and( isteps == 210 )):
        hadapp.applyAction(loadshedact)
        hadapp.applyAction(loadshedact2)
        hadapp.applyAction(loadshedact3)
        
    if (bloadshedding and( isteps == 215 )):
        hadapp.applyAction(loadshedact)
        hadapp.applyAction(loadshedact2)
        hadapp.applyAction(loadshedact3)
        
    if (bloadshedding and( isteps == 220 )):
        hadapp.applyAction(loadshedact)
        hadapp.applyAction(loadshedact2)
        hadapp.applyAction(loadshedact3)

    # execute one simulation time step	
    hadapp.executeDynSimuOneStep()       

    # check if it is the time to get the observations, 
    # outputob_nsimustep defines every outputob_nsimustep we get observation from the simulator
    if 	isteps%outputob_nsimustep == 0:
        before_getob_time = time.time()	
        ob_vals = hadapp.getObservations()
        after_getob_time = time.time()	
        total_dataconv_time += (after_getob_time - before_getob_time)
        
        '''  
        Implement under frequency load shedding here
        '''
    
        for itmp, ibus in enumerate(loadshedingbuslist):
            freqtmp = ob_vals[busfreqidx_dict[ibus]]
            loadorgvalue = loadorgvaluelist[itmp]
            
            #if the frequency of the bus is lower than the threshold
            if (freqtmp < freqthre and bUFLS and loadremainingvaluelist[itmp]>loadshedperc*0.5):
                if loadshedingtimer[itmp] == loadshedingtimersteps:  # if the load shedding delay timer is reached
                    
                    #define the load shedding action
                    loadshedact = gridpack.hadrec.Action()
                    loadshedact.actiontype = 3;
                    loadshedact.bus_number = ibus;
                    loadshedact.componentID = "1";
                    
                    # shed load amount should be the min value of load drop percentage and the remaining load value
                    loadshedamount = min(loadorgvalue*loadshedperc, loadremainingvaluelist[itmp])
                    loadshedact.percentage = -loadshedamount;  
                    hadapp.applyAction(loadshedact)
					
                    loadshedingtimer[itmp] = -1 # reset the timer once the load shedding action is applied
                    loadremainingvaluelist[itmp] -= loadshedamount # change the remaining load value as the load has been shed
                    print ('------Warning: Time: %5.3f seconds, bus %d frequency is %4.2f, shed %6.2f load at this bus, remaining load: %6.2f'%
                       (isteps*simu_time_step, ibus, freqtmp, loadshedamount, loadremainingvaluelist[itmp]))
                    
                else: # if the load shedding delay timer is not reached, increase the timer

                    loadshedingtimer[itmp] += 1
                
            else: # if the frequency of the bus is higher than the threshold, reset the timer
                loadshedingtimer[itmp] = -1
                            
        # compute line flows here
        lineflowonesteplist = []
        for iline in range(len(linedf)):
            brfrombus  =  linedf['From Bus  Number'][iline]# from bus number
            brtobus = linedf['To Bus  Number'][iline]
            cktid = str(linedf['Id'][iline])
            bus1mag = ob_vals[busmagidx_dict[brfrombus]] # get the line from bus voltage magnitude from the observation list
            bus2mag = ob_vals[busmagidx_dict[brtobus]] # get the line to bus voltage magnitude from the observation list
            bus1ang = ob_vals[busangidx_dict[brfrombus]] # get the line from bus voltage angle from the observation list, unit in rads
            bus2ang = ob_vals[busangidx_dict[brtobus]]  # get the line to bus voltage angle from the observation list, unit in rads
            bus1_volt_complx = complex (bus1mag*math.cos(bus1ang), bus1mag*math.sin(bus1ang))
            bus2_volt_complx = complex (bus2mag*math.cos(bus2ang), bus2mag*math.sin(bus2ang))
            
            y11 = Y11[iline]
            y12 = Y12[iline]
            y21 = Y21[iline]
            y22 = Y22[iline]
            
            # I=Y*V to compute the current of the line, for both from and to ends
            cur_from_complex = y11*bus1_volt_complx + y12*bus2_volt_complx 
            cur_to_complex = y21*bus1_volt_complx + y22*bus2_volt_complx 
            
            #print ('bus1_volt_complx: ', bus1_volt_complx)
            #print ('cur_from_complex: ', cur_from_complex)
            #print ('cur_from_complex.conjugate: ', cur_from_complex.conjugate())
            
            # Apparent power S = V*conj(I)
            S_from = bus1_volt_complx*cur_from_complex.conjugate()
            #print ('S_from: ', S_from)
            
            P_from = S_from.real
            Q_from = S_from.imag
            S_to = bus2_volt_complx*cur_to_complex.conjugate()
            P_to = S_from.real
            Q_to = S_from.imag
            
            #check whether the line flow S exceeds the line rating
            lineS = max(abs(S_from), abs(S_to))
            lineRating = lineparalist[iline][8]
            
            # if the line flow exceeds the rating, trip the line
            # lineS is with P.U. While Rating is actual value
            if ( blinetriprelay and lineS*100.0 > lineRating and nlinestatus[iline]==1 ): 
                
                if linetriptimer[iline] == linetriptimersteps:
                    #define a line tripping action here and apply it in the simulation
                    linetripact = gridpack.hadrec.Action()
                    linetripact.actiontype = 1;
                    linetripact.brch_from_bus_number = brfrombus;
                    linetripact.brch_to_bus_number = brtobus;
                    linetripact.componentID = cktid;
                
                    hadapp.applyAction(linetripact)
                
                    # set corresponding line status in the line status list to be 0, 
                    # if the line is tripped
                    nlinestatus[iline] = 0 
                
                    print ('------Warning: Time: %5.3f seconds, Branch from bus %d to bus %d with CKT %s has a flow %7.3f MVA, larger than the rating A %7.3f, line tripped!'%
                       (isteps*simu_time_step, brfrombus, brtobus, cktid, lineS*100.0, lineRating))
                else:
                    linetriptimer[iline] += 1
            else:
                linetriptimer[iline] = -1
            
            # if the line is tripped, line flow is zero                
            if nlinestatus[iline] == 0: 
                P_from = 0.0
                Q_from = 0.0
                
            # append the computed line flow in the list    
            lineflowonesteplist.append(P_from)
            lineflowonesteplist.append(Q_from)
        
        # add time in the head of the observation data list
        ob_vals.insert(0,isteps*simu_time_step)	
        observation_list.append(ob_vals)
        
        # add time in the head of the line flow data list
        lineflowonesteplist.insert(0,isteps*simu_time_step)
        lineflow_list.append(lineflowonesteplist)
        
        # add time in the head of the load data list
        loadvalue_onestep = loadremainingvaluelist[:]
        loadvalue_onestep.insert(0,isteps*simu_time_step)
        loadvalue_list.append(loadvalue_onestep)
        
    isteps = isteps + 1

gridpack_endtime = time.time()
print ("-------------!!! total gridpack time :  ", gridpack_endtime - gridpack_starttime)
print ("-------------!!! total gridpack powerflow and initialization time :  ", gridpack_initime - gridpack_starttime)
print ("-------------!!! total gridpack data converstion time :  ", total_dataconv_time)

# write the observations into a csv file
before_csvwr_time = time.time()
np_data = np.array(observation_list)

pd_data = pd.DataFrame(np_data,columns=csvhead)

csv_wrt_f = '14bus_test_example_observations'

pd_data.to_csv(csv_wrt_f+'.csv', index=False)

# also write line flow data into a csv file
np_lineflow_data = np.array(lineflow_list)

pd_lineflow_data = pd.DataFrame(np_lineflow_data,columns=lineflow_csvhead)

csv_lineflow_wrt_f = '14bus_test_example_lineflow'

pd_lineflow_data.to_csv(csv_lineflow_wrt_f+'.csv', index=False)

# also write load data into a csv file
np_load_data = np.array(loadvalue_list)

pd_load_data = pd.DataFrame(np_load_data,columns=loadvalue_csvhead)

csv_load_wrt_f = '14bus_test_example_load'

pd_load_data.to_csv(csv_load_wrt_f+'.csv', index=False)

total_csvwr_time = time.time() - before_csvwr_time
print("-------------!!! total csv write time :  ", total_csvwr_time)

# It's important to force the deallocation order here

hadapp = None
run_gridpack_necessary_env = None
