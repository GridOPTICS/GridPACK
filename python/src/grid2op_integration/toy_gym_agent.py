import numpy as np
from datetime import timedelta

import gymnasium
# Starting with 1.9.1 grid2op switch (as advised) from the legacy “gym” framework to the new “gymnasium” framework 
# (gym is no longer maintained since v0.26.2, see https://www.gymlibrary.dev/).

import sys
sys.path.append("/qfs/projects/gridpack_wind/grid2op_interface/grid2op_local/")
import grid2op
from grid2op.gym_compat import GymEnv

sys.path.append("/qfs/projects/gridpack_wind/grid2op_interface/GridPACK/python/src/")
from grid2op_backend import GridPACKBackend
NUMBA_ = False

# create grid2op environment
env = grid2op.make(
    dataset="/qfs/projects/gridpack_wind/grid2op_interface/GridPACK/python/src/grid2op_integration/test_grid2op",
    grid_path="input_9bus.xml",
    backend=GridPACKBackend( # GridPACK Backend for Grid2Op
        log_freq=1, # at what frequency to log the data
        gridpack_stepsize=0.005, # step size for GridPACK
        grid2op_stepsize=5, # step size for grid2op
        detailed_infos_for_cascading_failures=False,
        can_be_copied=True,
        grid_path="input_9bus.xml"
    ),
    data_feeding_kwargs={
        # start datatime goes here
        # for frequency - a pull request from Grid2Op
        "time_interval": timedelta(seconds=5)
    },
    # experimental_read_from_local_dir=True
)
# env = grid2op.make("l2rpn_case14_sandbox")
# reset the environment
print("[INFO] Grid2Op environment created, will reset the environment now.")
obs = env.reset()
print("[INFO] Environment reset, will now convert this into a Gym environment.")

# create the gymnasium environment
gym_env = GymEnv(env)  
print(f"The \"gym_env\" is a gym environment: {isinstance(gym_env, gymnasium.Env)}")

dim_act_space = np.sum([np.sum(gym_env.action_space[el].shape) for el in gym_env.action_space.spaces])
print(f"The size of the action space is : "
      f"{dim_act_space}")
dim_obs_space = np.sum([np.sum(gym_env.observation_space[el].shape).astype(int) 
                        for el in gym_env.observation_space.spaces])
print(f"The size of the observation space is : "
      f"{dim_obs_space}")
# breakpoint()
obs_gym, info = gym_env.reset()

# do nothing
gym_act = {}
obs, reward, done, truncated, info = gym_env.step(gym_act)
print(f"Reward: {reward}")

gym_env.close()
# env.close()
