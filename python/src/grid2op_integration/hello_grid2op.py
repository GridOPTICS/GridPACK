import warnings
import sys, time
from memory_profiler import memory_usage
from datetime import timedelta

sys.path.append("/qfs/projects/gridpack_wind/grid2op_interface/grid2op_local/")
import grid2op
from grid2op.PlotGrid.PlotMatplot import PlotMatplot
from grid2op.Chronics import ChangeNothing

from agent import LoadSheddingAgent, DoNothingAgent
from utils import parse_arguments, DataCollector

sys.path.append("/qfs/projects/gridpack_wind/grid2op_interface/GridPACK/python/src/")
from grid2op_backend import GridPACKBackend

def main():
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
                            gridpack_stepsize=0.005,
                            grid2op_stepsize=args.grid2op_stepsize,
                            can_be_copied=True,
                            grid_path=filename
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
    # dc.log_data(env)
    obs = env.get_obs()
    print(obs.load_p)
    breakpoint()

    # reset environment
    print("============ Environment Reset =============")
    obs = env.reset()
    print(obs.load_p)
    # log data
    dc.log_data(env)
    
    reward = env.reward_range[0]
    done = False

    if args.shed_load:
        my_agent = LoadSheddingAgent(env.action_space)
    else:
        my_agent = DoNothingAgent(env.action_space)
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
            if args.shed_load:
                filename_prefix += "_lsagent"
            else:
                filename_prefix += "_donothingagent"
            dc.save_data(filename_prefix)

        # update counter
        counter += 1

    # print(env.observation_space)
    print(f"[INFO] Runtime: {time.perf_counter()-start_time} seconds")

    # plot_helper = PlotMatplot(env.observation_space)
    # fig = plot_helper.plot_layout()

if __name__=="__main__":
    # Measure memory and time
    start_time = time.perf_counter()
    mem_usage = memory_usage((main, ), max_usage=True)
    end_time = time.perf_counter()

    print(f"Total runtime: {end_time - start_time:.4f} seconds")
    print(f"Peak memory usage: {mem_usage:.2f} MiB")
    # main()