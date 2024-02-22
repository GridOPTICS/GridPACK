import numpy as np
from datetime import timedelta

import gymnasium
# Starting with 1.9.1 grid2op switch (as advised) from the legacy “gym” framework to the new “gymnasium” framework 
# (gym is no longer maintained since v0.26.2, see https://www.gymlibrary.dev/).

import grid2op
from grid2op.gym_compat import GymEnv

from agent import LoadSheddingAgent
from reward import HADRECReward

import sys
sys.path.append("/qfs/projects/gridpack_wind/grid2op_interface/GridPACK/python/src/")
from grid2op_backend import GridPACKBackend
NUMBA_ = False

hadrec_reward = HADRECReward(preFaultActionPenalty=1, actionPenalty=2, invalidActionPenalty=42)

# create grid2op environment
env = grid2op.make(
    dataset="/qfs/projects/gridpack_wind/grid2op_interface/GridPACK/python/src/grid2op_integration/test_grid2op",
    grid_path="input_9bus.xml",
    backend=GridPACKBackend( # GridPACK Backend for Grid2Op
        log_freq=1, # at what frequency to log the data
        grid2op_stepsize=5, # step size for grid2op
        detailed_infos_for_cascading_failures=False,
        can_be_copied=True
    ),
    data_feeding_kwargs={
        # start datatime goes here
        # for frequency - a pull request from Grid2Op
        "time_interval": timedelta(seconds=5)
    },
    reward_class=hadrec_reward,
    # experimental_read_from_local_dir=True
)

# reset the environment
print("[INFO] Grid2Op environment created, will reset the environment now.")
obs = env.reset()
print("[INFO] Environment reset, will now convert this into a Gym environment.")

# create the gymnasium environment
gym_env = GymEnv(env)  
print(f"The \"gym_env\" is a gym environment: {isinstance(gym_env, gymnasium.Env)}")

# Reset the gymnasium environment
obs_gym, info = gym_env.reset()

# do nothing
gym_act = {}
obs, reward, done, truncated, info = gym_env.step(gym_act)
print(f"Reward: {reward}")

# gym_env.close()
# here we run one episode without any load shedding actions
obs_noact = list()
actions = list()
episode_rew = 0.0
rollout_length = 10

my_agent = LoadSheddingAgent(env.action_space)
for cnt in range(rollout_length):
    obs_noact.append(obs)
    
    # simulate the env by taking the actions from action_lst to the next step
    # and get the observations and reward for this step
    # action = my_agent.act(obs, rew, done)
    obs, rew, done, _, _ = gym_env.step(gym_act)
    print(obs["gen_p"])

    episode_rew += rew

    # check whether this episode is done
    if done:
        break

print("------------------- Total steps: %d, Episode total reward without any load shedding actions: "%(cnt), episode_rew)

# How to capture observations? - access it like a dictionary - this page lists all the information that you can capture from observation - https://grid2op.readthedocs.io/en/latest/user/observation.html
# How to update actions? - we can do that using grid2op agent and the action gets executed in the backend. 
# How to define reward function? - you can define your own reward class
# How to define done?
# How to train policy?