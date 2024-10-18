import grid2op
import warnings
import sys
sys.path.append("/qfs/projects/gridpack_wind/grid2op_interface/GridPACK/python/src/")
from grid2op_backend import GridPACKBackend
from grid2op.PlotGrid.PlotMatplot import PlotMatplot
from grid2op.Agent import RandomAgent

import argparse

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--gridpack_config", 
        help="Configuration file to run the simulation.", 
        default="input_9b3g.xml"
    )
    parser.add_argument(
        "--grid2op_config", 
        help="Configuration file for Grid2Op.", 
        default="/qfs/projects/gridpack_wind/grid2op_interface/GridPACK/python/src/example/test_grid2op"
    )
    parser.add_argument(
        "--iter", 
        help="Total number of iterations to run the code.", 
        default=200,
        type=int
    )
    parser.add_argument(
        "--log_freq", 
        help="Number of steps between each print of the data.", 
        default=1,
        type=int
    )

    return parser.parse_args()
    
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
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        env = grid2op.make(config_filepath,
                        grid_path=filename,
                        backend=GridPACKBackend(save_at=args.iter-1, log_freq=args.log_freq)
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
    while counter < args.iter:
        new_load_p = obs.load_p * 1.1
        new_load_q = obs.load_q * 0.9
        # print(obs.load_p, obs.load_q)
        action = env.action_space(
            {"injection": {
                "load_p": new_load_p,
                "load_q": new_load_q
                }
            }
        )
        obs, reward, done, info = env.step(action)
        counter += 1

    print(env.observation_space)

    # plot_helper = PlotMatplot(env.observation_space)
    # fig = plot_helper.plot_layout()
