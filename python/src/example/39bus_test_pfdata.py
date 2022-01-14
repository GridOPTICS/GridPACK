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

print ('before gridpack ini')
noprintflag = gridpack.NoPrint()
noprintflag.setStatus(True) # if you set this to be True, no any print out from GridPACK Hadrec module, usefull for parallel running with ray to disable any temp print outs
run_gridpack_necessary_env = gridpack.Environment()
hadapp = gridpack.hadrec.Module()

print ('after gridpack ini')

gridpack_starttime = time.time()


hadapp.readPowerFlowData(inname, 0)  

print ('------first read load par')
loadbus = 508
loadid = '1'
parstr = 'LOAD_PL'
testpar = hadapp.getDataCollectionLoadParam(loadbus, loadid, parstr)
print('LOAD_PL:', testpar)

parstr = 'LOAD_QL'
testpar = hadapp.getDataCollectionLoadParam(loadbus, loadid, parstr)
print('LOAD_QL:', testpar)

# modify the P and Q of a specific load
print ('\n----------test again, first modify the load \n')
loadbus = 508
loadid = '1'
loadpar = 'LOAD_PL'
newloadp = 600.0  # this is the new value of the load P, unit is MW
hadapp.modifyDataCollectionLoadParam(loadbus, loadid, loadpar, newloadp)

loadpar = 'LOAD_QL'
newloadq = 300.0  # this is the new value of the load Q, unit is MVar
hadapp.modifyDataCollectionLoadParam(loadbus, loadid, loadpar, newloadq)

parstr = 'LOAD_PL'
testpar = hadapp.getDataCollectionLoadParam(loadbus, loadid, parstr)
print('LOAD_PL:', testpar)

parstr = 'LOAD_QL'
testpar = hadapp.getDataCollectionLoadParam(loadbus, loadid, parstr)
print('LOAD_QL:', testpar)

genbus = 38
genid = '1'
genpar = 'GENERATOR_PG'
testpar = hadapp.getDataCollectionGenParam(genbus, genid, genpar)
print('GENERATOR_PG:', testpar)

genpar = 'GENERATOR_PF_PGEN'
testpar = hadapp.getDataCollectionGenParam(genbus, genid, genpar)
print('GENERATOR_PF_PGEN:', testpar)

genpar = 'GENERATOR_MBASE'
testpar = hadapp.getDataCollectionGenParam(genbus, genid, genpar)
print('GENERATOR_MBASE:', testpar)

print ('now solve power flow ')
# solve the power flow
hadapp.solvePowerFlow()   

print ('after solve power flow ')

genbus = 38
genid = '1'
genpar = 'GENERATOR_PG'
testpar = hadapp.getDataCollectionGenParam(genbus, genid, genpar)
print('GENERATOR_PG:', testpar)

genpar = 'GENERATOR_PF_PGEN'
testpar = hadapp.getDataCollectionGenParam(genbus, genid, genpar)
print('GENERATOR_PF_PGEN:', testpar)                                                                                                                                                                











