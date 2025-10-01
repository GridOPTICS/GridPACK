import sys
sys.path.append("../")
import numpy as np
import pandas as pd

import grid2op_backend

# initiate backend
backend = grid2op_backend.GridPACKBackend()
# load grid
backend.load_grid("input_39bus_step005_v33.xml")

# observation counter
count = 1

# get the observation name list 
(obs_genBus, obs_genIDs, obs_loadBuses, obs_loadIDs, obs_busIDs) = backend._hadapp.getObservationLists()

# create observation names for csv file header writting purpose
csvhead = []
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

# get observations
ob_vals = backend._hadapp.getObservations()
print("Observation (%d, %s):" %(count, ob_vals))

# store observations
observation_list = []

# while simulation is working
while (not backend._hadapp.isDynSimuDone()):
    # run simulation
    backend.runpf()

    # get observation
    ob_vals = backend._hadapp.getObservations()
    print("Observation (%d, %s):" %(count, ob_vals))
    
    # append to list 
    observation_list.append(ob_vals)
    
    # increment the counter
    count += 1

# store to numpy array
np_data = np.array(observation_list)
# convert to pd.DataFrame
pd_data = pd.DataFrame(np_data,columns=csvhead)

# csv file
csv_wrt_f = '39bus_test_grid2op_observations'
# store csv
pd_data.to_csv(csv_wrt_f+'.csv', index=False)

backend.close()