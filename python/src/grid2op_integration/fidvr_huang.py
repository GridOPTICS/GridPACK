import sys
import numpy as np
from datetime import timedelta

sys.path.append("/qfs/projects/gridpack_wind/grid2op_interface/grid2op_local/")
import grid2op
from grid2op import make
from grid2op.Action import BaseAction
from grid2op.Reward import BaseReward
from grid2op.gym_compat import GymEnv, BoxGymObsSpace, DiscreteActSpaceGymnasium, DiscreteActSpace

import stable_baselines3
from stable_baselines3 import DQN
from stable_baselines3.common.vec_env import DummyVecEnv, VecFrameStack

import gymnasium as gym

# Print versions
print(f"Grid2Op Version: {grid2op.__version__}")
print(f"Stable Baseline3 Version: {stable_baselines3.__version__}")
print(f"Gymnasium Version: {gym.__version__}")

sys.path.append("/qfs/projects/gridpack_wind/grid2op_interface/GridPACK/python/src/")
from grid2op_backend import GridPACKBackend
from test_probe import test_probe_run_gym

# Constants
TVRC_THRESHOLDS = [
    (0.33, 0.7),
    (0.5, 0.8),
    (1.5, 0.9),
    (np.inf, 0.95)
]
T_FAULT_CLEAR = 2.0  # seconds (assumed time of fault clearance)
T_FINAL_CHECK = T_FAULT_CLEAR + 4.0  # 4 seconds after fault clearance
LOAD_SHED_RATIO = 0.2  # 20% load shedding
C1, C2, C3 = 1.0, 10.0, 100.0  # Reward weights

FILENAME = "input_9bus.xml"
CONFIG_FILEPATH = "/qfs/projects/gridpack_wind/grid2op_interface/GridPACK/python/src/grid2op_integration/test_grid2op"

# Controlled load buses (by load ID/index, not busbar)
# You may need to adjust these indices based on your grid
CONTROLLED_LOADS = [0, 1, 2]  # Example: first 3 loads; verify in your case


# === FIDVR Reward Class ===
class FIDVRReward(BaseReward):
    def __init__(self):
        super().__init__()
        self.init_loads = None
        self.previous_loads = None

    def __call__(self, action, env, has_error, is_done, is_illegal, is_ambiguous) -> float:
        obs = env.get_obs()
        if self.init_loads is None and obs is not None:
            self.init_loads = obs.load_p.copy()
            self.previous_loads = self.init_loads.copy()

        voltages = obs.load_v
        loads = obs.load_p
        time = obs.current_step * env._time_step  # seconds

        if time == 0.0:
            self.previous_loads = self.init_loads.copy()
            return 0.0

        if time > T_FINAL_CHECK and (voltages < 0.95).any():
            return -1000.0

        v_penalty = 0.0
        if time > T_FAULT_CLEAR:
            t_elapsed = time - T_FAULT_CLEAR
            if t_elapsed < 0.33:
                v_penalty += np.sum(np.minimum(voltages - 0.7, 0.0))
            elif t_elapsed < 0.5:
                v_penalty += np.sum(np.minimum(voltages - 0.8, 0.0))
            elif t_elapsed < 1.5:
                v_penalty += np.sum(np.minimum(voltages - 0.9, 0.0))
            else:
                v_penalty += np.sum(np.minimum(voltages - 0.95, 0.0))

        # p.u. penalty
        eps = 1e-6
        load_diff = (self.previous_loads - loads) / np.maximum(self.init_loads, eps)
        load_penalty = np.sum(np.maximum(load_diff, 0.0))
        self.previous_loads = loads.copy()

        # 1 per invalid bus
        invalid_shedding = (loads < 1e-3) & (load_diff > 1e-6)
        invalid_penalty = float(np.sum(invalid_shedding))

        return float(C1 * v_penalty - C2 * load_penalty - C3 * invalid_penalty)


# === Create Environment ===
env = make(
    CONFIG_FILEPATH,
    grid_path=FILENAME,
    backend=GridPACKBackend(
        log_freq=1,
        grid2op_stepsize=0.05,  # 50ms internal steps
        can_be_copied=True
    ),
    data_feeding_kwargs={
        "time_interval": timedelta(seconds=0.05)  # 50ms per step
    },
    reward_class=FIDVRReward,
    action_class=BaseAction
)

# Set up Gym environment
gym_env = GymEnv(env)

# Define observation space: only voltage and load_p
gym_env.observation_space = BoxGymObsSpace(
    env.observation_space,
    attr_to_keep=["load_v", "load_p"]
)

# === Define Discrete Action Space: 2^n combinations ===
n_controlled = len(CONTROLLED_LOADS)
init_loads = env.get_obs().load_p

STAGES = [0.0, 0.2, 0.4]   # 3-stage UVLS per bus (0%, 20%, 40% of P0)
# total actions = (len(STAGES)) ** n_controlled

# === Define Discrete Action Space: 2^n combinations via `injection` ===
n_controlled = len(CONTROLLED_LOADS)
init_obs = env.get_obs()
init_p = init_obs.load_p.copy()
init_q = init_obs.load_q.copy()

action_list = []
for mask in range(2 ** n_controlled):
    new_p = init_p.copy()
    new_q = init_q.copy()
    for bit, load_id in enumerate(CONTROLLED_LOADS):
        if (mask >> bit) & 1:
            frac = 1.0 - LOAD_SHED_RATIO          # 0.8 for 20% shed
            new_p[load_id] = frac * init_p[load_id]
            new_q[load_id] = frac * init_q[load_id]  # preserve power factor
    act = env.action_space({
        "injection": {
            "load_p": new_p,
            "load_q": new_q
        }
    })
    action_list.append(act)

# Bind the discrete action space (keep the class you already use)
gym_env.action_space = DiscreteActSpace(env.action_space, action_list=action_list)

# -- Sanity poke: does ALL_ONES action change load_p? --
obs0 = gym_env.reset()
ALL_ONES_IDX = (1 << len(CONTROLLED_LOADS)) - 1

# get the last-frame vector (handles VecFrameStack flattening)
obs_vec = obs0[0]
per = len(obs_vec) // getattr(gym_env, "n_stack", 1) if hasattr(gym_env, "n_stack") else len(obs_vec)
frame0 = obs_vec[-per:]
half = per // 2
print([ALL_ONES_IDX] if hasattr(gym_env,"num_envs") else ALL_ONES_IDX)
# because attr_to_keep=["v","load_p"] now:
load_p0 = frame0[half:]#.copy()

obs1, reward1, done1, x, info1 = gym_env.step([ALL_ONES_IDX] if hasattr(gym_env,"num_envs") else ALL_ONES_IDX)

obs_vec1 = obs1[0] if hasattr(obs1, "ndim") and obs1.ndim > 1 else obs1
frame1 = obs_vec1[-per:]
load_p1 = frame1[half:]#.copy()

print("load_p before:", load_p0)
print("load_p after :", load_p1)

# sys.exit(1)

# After you’ve built gym_env (and set the discrete action space)
# After gym_env is ready (and after you set the discrete action space)
test_probe_run_gym(
    gym_env,
    CONTROLLED_LOADS,
    max_steps=60,
    T_FAULT_CLEAR=2.0,
    ts=0.05,
    n_stack=4,
    obs_order="vp"   # per-frame is [load_v, load_p]; switch to "pv" if your numbers look like [load_p, load_v]
)


# sys.exit(1)

# === Stack Observations (History) ===
# Required: VecEnv wrapper for VecFrameStack
gym_env = DummyVecEnv([lambda: gym_env])
gym_env = VecFrameStack(gym_env, n_stack=4)  # Use last 4 observations

# === Initialize DQN Model ===
model = DQN(
    "MlpPolicy",
    gym_env,
    learning_rate=1e-3,
    buffer_size=10000,
    learning_starts=200,
    batch_size=32,
    train_freq=4,
    target_update_interval=100,
    verbose=1,
    tensorboard_log="./dqn_fidvr_tensorboard/"
)

# === Train the Model ===
print("Starting training...")
model.learn(total_timesteps=5, log_interval=1)

# Save model
model.save("dqn_fidvr_load_shedding")

print("Training complete. Model saved.")