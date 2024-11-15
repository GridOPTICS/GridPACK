import grid2op
import warnings
import sys, time
from datetime import timedelta
import pandas as pd

sys.path.append("/qfs/projects/gridpack_wind/grid2op_interface/GridPACK/python/src/")
from grid2op_backend import GridPACKBackend
from grid2op.PlotGrid.PlotMatplot import PlotMatplot
from grid2op.Agent import RandomAgent
from grid2op.Chronics import ChangeNothing

import argparse

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--gridpack_config", 
        help="Configuration file to run the simulation.", 
        default="input_9bus.xml",
        choices=["input_9bus.xml", "input_39bus_IBR.xml"]
    )
    parser.add_argument(
        "--grid2op_config", 
        help="Configuration file for Grid2Op.", 
        default="/qfs/projects/gridpack_wind/grid2op_interface/GridPACK/python/src/example/test_grid2op"
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
    
BUS_LOGGER = []
GEN_LOGGER = []
LOAD_LOGGER = []
BRANCH_LOGGER = []

def save_data(filename_prefix):
    print("[INFO] Saving the data")
    # bus data
    res_bus = pd.concat(BUS_LOGGER, ignore_index=True)
    # adjusting for counter increase during warm up calls
    res_bus["tick"] = res_bus["tick"] - res_bus["tick"].values[0]
    res_bus.to_csv(f"/qfs/projects/gridpack_wind/grid2op_interface/temp/{filename_prefix}_res_bus.csv")
    
    # gen data
    res_gen = pd.concat(GEN_LOGGER, ignore_index=True)
    # adjusting for counter increase during warm up calls
    res_gen["tick"] = res_gen["tick"] - res_gen["tick"].values[0]
    res_gen.to_csv(f"/qfs/projects/gridpack_wind/grid2op_interface/temp/{filename_prefix}_res_gen.csv")
    
    # load data
    res_load = pd.concat(LOAD_LOGGER, ignore_index=True)
    # adjusting for counter increase during warm up calls
    res_load["tick"] = res_load["tick"] - res_load["tick"].values[0]
    res_load.to_csv(f"/qfs/projects/gridpack_wind/grid2op_interface/temp/{filename_prefix}_res_load.csv")
    
    # line data
    res_line = pd.concat(BRANCH_LOGGER, ignore_index=True)
    # adjusting for counter increase during warm up calls
    res_line["tick"] = res_line["tick"] - res_line["tick"].values[0]
    res_line.to_csv(f"/qfs/projects/gridpack_wind/grid2op_interface/temp/{filename_prefix}_res_line.csv")

if __name__=="__main__":
    # parse arguments
    args = parse_arguments()
    
    # NOTE: You must run this code from the application folder
    # filename = "input_240bus.xml"
    # filename = "input_9b3g.xml"
    # filename = "input_39bus_step005_v33.xml" ## formatting issues
    filename = args.gridpack_config
    
    # config_filepath = "/qfs/projects/gridpack_wind/grid2op_interface/GridPACK/python/src/example/test_grid2op"
    config_filepath = args.grid2op_config

    # make environment
    print("============ Initializing Environment =============")
    start_time = time.perf_counter()
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        env = grid2op.make(config_filepath,
                        grid_path=filename,
                        backend=GridPACKBackend(
                            log_freq=args.log_freq,
                            grid2op_stepsize=args.grid2op_stepsize
                        ),
                        data_feeding_kwargs={
                            "time_interval": timedelta(
                                seconds=args.grid2op_stepsize
                            )}
                    ) # mention the time resolution - resolution of the environment - episode size goes - stuff from config.py file won't be used if mentioned here - episode can be updated in the config files

    # reset environment
    print("============ Environment Reset =============")
    # agent = RandomAgent(env.action_space)
    obs = env.reset()
        
    # print network analytics
    print("Number of buses:  %d" % (obs.n_sub))
    print("Number of generators: %d" % (obs.n_gen))
    print("Number of loads: %d" % (obs.n_load))
    print("Number of lines: %d" % (obs.n_line))
    print("Number of storage units: %d" % (obs.n_storage))

    # step environment
    counter = 0
    while counter < args.grid2op_steps:
        new_load_p = obs.load_p * 1.1
        new_load_q = obs.load_q * 1.1
        # print(obs.load_p, obs.load_q)
        action = env.action_space(
            {
                "injection": {
                    "load_p": new_load_p,
                    "load_q": new_load_q
                }
            }
        )
        # action = env.action_space({})

        # start time to save the data finally
        grid2op_start_time = env.time_stamp
        # print(env.time_stamp, env.backend._counter)
        # print(env.delta_time_seconds)
        obs, reward, done, info = env.step(action)

        # log data
        BUS_LOGGER += env.backend.bus_logger
        GEN_LOGGER += env.backend.gen_logger
        LOAD_LOGGER += env.backend.load_logger
        BRANCH_LOGGER += env.backend.branch_logger
        
        # save data
        if counter == args.grid2op_steps-1: # save at 2nd grid2op step
            filename_prefix = filename.split(".")[0]
            save_data(filename_prefix)

        # update counter
        counter += 1

    # print(env.observation_space)
    print(f"[INFO] Runtime: {time.perf_counter()-start_time} seconds")

    # plot_helper = PlotMatplot(env.observation_space)
    # fig = plot_helper.plot_layout()
