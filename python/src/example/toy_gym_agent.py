import numpy as np
from datetime import timedelta

import grid2op
from grid2op.gym_compat import GymEnv

import sys
sys.path.append("/qfs/projects/gridpack_wind/grid2op_interface/GridPACK/python/src/")
from grid2op_backend import GridPACKBackend
NUMBA_ = False

# create grid2op environment
env = grid2op.make(
    dataset="/qfs/projects/gridpack_wind/grid2op_interface/GridPACK/python/src/example/test_grid2op",
    grid_path="input_9bus.xml",
    backend=GridPACKBackend( # GridPACK Backend for Grid2Op
        log_freq=1, # at what frequency to log the data
        grid2op_stepsize=5 # step size for grid2op
    ),
    data_feeding_kwargs={
        # start datatime goes here
        # for frequency - a pull request from Grid2Op
        "time_interval": timedelta(seconds=5)
    }
)
breakpoint()
# env = grid2op.make("l2rpn_case14_sandbox")

# Starting with 1.9.1 grid2op switch (as advised) from the legacy “gym” framework to the new “gymnasium” framework 
# (gym is no longer maintained since v0.26.2, see https://www.gymlibrary.dev/).

# print(env.observation_space.gen_p)
# env.observation_space["gen_p"].low[:] = -np.inf
# env.observation_space["gen_p"].high[:] = np.inf

gym_env = GymEnv(env)  # create the gymnasium environment
# obs = env.reset()