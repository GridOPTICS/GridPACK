# ------------ Config ------------
# GRID_XML = "input_9bus.xml"
# GRID_XML = "input_39bus_IBR.xml"
# GRID_XML = "input_145.xml"

#### failing
# GRID_XML = "input_240.xml"
# GRID_XML = "input_3000.xml"

#### different formatting
# GRID_XML = "input_39bus_step005_v33.xml"

CONFIG_DIR = f"/qfs/projects/gridpack_wind/grid2op_interface/GridPACK/python/src/grid2op_integration"
OBS_ATTRS = ["load_v", "load_p"]  # used for the NN state; logging handles more attrs when present
CONTROLLED_LOADS = {
    "input_9bus.xml": [4, 5, 7, 9],
    "input_9bus_fault.xml": [4, 5, 7, 9],
    # "input_39bus_step005_v33.xml": [38, 39],
    "input_39bus_IBR.xml": [2, 3, 6, 7, 11, 14, 15, 17, 19, 20, 22, 23, 24, 25, 26, 27, 28, 30, 38],
    "input_39bus_IBR_fault.xml": [2, 3, 6, 7, 11, 14, 15, 17, 19, 20, 22, 23, 24, 25, 26, 27, 28, 30, 38],
    "input_240.xml": [31],
    "input_145.xml": [33, 34, 50, 57, 65, 67, 69, 70, 73, 77, 78, 79, 80, 81, 83, 84, 87, 88, 89, 91, 92, 93, 94, 98, 100, 101, 103, 104, 105, 106, 109, 110, 111, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144],
    "input_145_fault.xml": [33, 34, 50, 57, 65, 67, 69, 70, 73, 77, 78, 79, 80, 81, 83, 84, 87, 88, 89, 91, 92, 93, 94, 98, 100, 101, 103, 104, 105, 106, 109, 110, 111, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144]
}
SHED_RATIO = 0.2
TIME_STEP_S = 1
TOTAL_STEPS = 30000
SAVE_EVERY_STEPS = 10
TARGET_SYNC = 1000
BATCH_SIZE = 4
GAMMA = 0.99
LR = 1e-3
EPS_START, EPS_END, EPS_DECAY = 1.0, 0.05, 15000

# ------------ Reward ------------
T_FAULT_CLEAR = 2.0
T_FINAL_CHECK = T_FAULT_CLEAR + 4.0
C1, C2, C3 = 260.0, 150.0, 3.0