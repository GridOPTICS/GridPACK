import argparse
import pandas as pd

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--gridpack_config", 
        help="Configuration file to run the simulation.", 
        default="input_9bus.xml",
        # choices=["input_9bus.xml", "input_39bus_IBR.xml", "input_240bus.xml", "input_240.xml"]
    )
    parser.add_argument(
        "--grid2op_config", 
        help="Configuration file for Grid2Op.", 
        default="/qfs/projects/gridpack_wind/grid2op_interface/GridPACK/python/src/grid2op_integration/test_grid2op"
    )
    parser.add_argument(
        "--shed_load", 
        help="Whether to deploy load shedding agent.", 
        action='store_true',
    )
    parser.add_argument(
        "--grid2op_stepsize", 
        help="Timestep size for Grid2Op (in seconds).", 
        default=1, # 1 minute NOTE: Update it to 60 later
        type=int
    )
    parser.add_argument(
        "--grid2op_steps", 
        help="Total number of iterations to run the code.", 
        default=5,
        type=int
    )
    parser.add_argument(
        "--log_freq", 
        help="Number of steps between each print of the data.", 
        default=1,
        type=int
    )

    return parser.parse_args()

class DataCollector():
    def __init__(self):
        self.BUS_LOGGER = []
        self.GEN_LOGGER = []
        self.LOAD_LOGGER = []
        self.LINE_LOGGER = []
        self.TRAFO_LOGGER = []

    def log_data(self, env):
        # log data
        self.BUS_LOGGER += env.backend.bus_logger
        self.GEN_LOGGER += env.backend.gen_logger
        self.LOAD_LOGGER += env.backend.load_logger
        self.LINE_LOGGER += env.backend.line_logger
        self.TRAFO_LOGGER += env.backend.trafo_logger

    def save_data(self, filename_prefix):
        print("[INFO] Saving the data")
        # bus data
        res_bus = pd.concat(self.BUS_LOGGER, ignore_index=True)
        print(res_bus.head())
        # adjusting for counter increase during warm up calls
        # res_bus["tick"] = res_bus["tick"] - res_bus["tick"].values[0]
        res_bus.to_csv(f"/qfs/projects/gridpack_wind/grid2op_interface/temp/{filename_prefix}_res_bus.csv")
        
        # gen data
        res_gen = pd.concat(self.GEN_LOGGER, ignore_index=True)
        # adjusting for counter increase during warm up calls
        # res_gen["tick"] = res_gen["tick"] - res_gen["tick"].values[0]
        res_gen.to_csv(f"/qfs/projects/gridpack_wind/grid2op_interface/temp/{filename_prefix}_res_gen.csv")
        
        # load data
        res_load = pd.concat(self.LOAD_LOGGER, ignore_index=True)
        # adjusting for counter increase during warm up calls
        # res_load["tick"] = res_load["tick"] - res_load["tick"].values[0]
        res_load.to_csv(f"/qfs/projects/gridpack_wind/grid2op_interface/temp/{filename_prefix}_res_load.csv")
        
        # line data
        res_line = pd.concat(self.LINE_LOGGER, ignore_index=True)
        # adjusting for counter increase during warm up calls
        # res_line["tick"] = res_line["tick"] - res_line["tick"].values[0]
        res_line.to_csv(f"/qfs/projects/gridpack_wind/grid2op_interface/temp/{filename_prefix}_res_line.csv")

        # trafo data
        res_trafo = pd.concat(self.TRAFO_LOGGER, ignore_index=True)
        # adjusting for counter increase during warm up calls
        # res_trafo["tick"] = res_trafo["tick"] - res_trafo["tick"].values[0]
        res_trafo.to_csv(f"/qfs/projects/gridpack_wind/grid2op_interface/temp/{filename_prefix}_res_trafo.csv")