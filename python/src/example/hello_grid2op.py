import warnings
import sys, time
import argparse

from datetime import timedelta
import pandas as pd
import pdb; 
# sys.path.append("/qfs/projects/gridpack_wind/grid2op_interface/install_local/grid2op")
import grid2op
from grid2op.PlotGrid.PlotMatplot import PlotMatplot
from grid2op.Agent import BaseAgent, RandomAgent
from grid2op.Chronics import ChangeNothing

sys.path.append("/qfs/projects/gridpack_wind/grid2op_interface/GridPACK/python/src/")
from grid2op_backend import GridPACKBackend

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--gridpack_config", 
        help="Configuration file to run the simulation.", 
        default="input_9bus.xml",
        choices=["input_9bus.xml", "input_39bus_IBR.xml", "input_240bus.xml", "input_240.xml"]
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


class LoadSheddingAgent(BaseAgent):
    def __init__(self, action_space):
        # define here the constructor of your agent
        # here we say our agent needs "something_else" and "and_another_something"
        # to be built just to demonstrate it does not cause any problem to extend the
        # construction of the base class BaseAgent that only takes "action_space" as a constructor
        BaseAgent.__init__(self, action_space)
        self.do_nothing = self.action_space({})

    def act(self, obs, reward, done=False):
        if ((obs.current_step >= 2) & (obs.current_step < 4)):
            new_load_p = obs.load_p * 1.2
            new_load_q = obs.load_q * 1.2
            
            # this is the only method you need to implement
            # it takes an observation obs (and a reward and a flag)
            # and should return a valid action
            dictionary_describing_the_action = {
                    "injection": {
                        "load_p": new_load_p,
                        "load_q": new_load_q
                    }
                }  # this can be anything you want that grid2op understands
            
            my_action = self.action_space(dictionary_describing_the_action)
        else:
            my_action = self.do_nothing
        return my_action
    

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
                            grid2op_stepsize=args.grid2op_stepsize,
                            can_be_copied=False
                        ),
                        data_feeding_kwargs={
                            # start datatime goes here
                            # for frequency - a pull request from Grid2Op
                            "time_interval": timedelta(
                                seconds=args.grid2op_stepsize
                            )}
                    ) # mention the time resolution - resolution of the environment - episode size goes - stuff from config.py file won't be used if mentioned here - episode can be updated in the config files
        # z-environment is passed to the backend - this will be available in the next version.
    
    dc = DataCollector()
    dc.log_data(env)
    
    # reset environment
    print("============ Environment Reset =============")
    # pdb.set_trace()
    obs = env.reset()
    # log data
    dc.log_data(env)
    
    reward = env.reward_range[0]
    done = False

    my_agent = LoadSheddingAgent(env.action_space)
    # print("======== Helloooooo ")
        
    # print network analytics
    print("Number of buses:  %d" % (obs.n_sub))
    print("Number of generators: %d" % (obs.n_gen))
    print("Number of loads: %d" % (obs.n_load))
    print("Number of lines: %d" % (obs.n_line))
    print("Number of storage units: %d" % (obs.n_storage))
    # print(obs.nb_substation)

    # step environment
    counter = 0
    while counter < args.grid2op_steps:
        # start time to save the data finally
        grid2op_start_time = env.time_stamp
        # need to mention grid2op start time - probably in chronics. start_datetime
        # print(env.time_stamp, env.backend._counter)
        # print(env.delta_time_seconds)
        action = my_agent.act(obs, reward, done)
        obs, reward, done, info = env.step(action)
        print(f"Reward: {reward}, Done: {done}")

        # log data
        dc.log_data(env)
        
        # save data
        if counter == args.grid2op_steps-1: # save at 2nd grid2op step
            filename_prefix = filename.split(".")[0]
            dc.save_data(filename_prefix)

        # update counter
        counter += 1

    # print(env.observation_space)
    print(f"[INFO] Runtime: {time.perf_counter()-start_time} seconds")

    # plot_helper = PlotMatplot(env.observation_space)
    # fig = plot_helper.plot_layout()