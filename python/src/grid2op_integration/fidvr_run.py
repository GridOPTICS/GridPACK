#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from bs4 import BeautifulSoup
import os, sys, json, time, math, uuid, argparse
from dataclasses import dataclass, asdict
from typing import List, Dict, Any, Optional
import numpy as np
from datetime import timedelta
import datetime as dt

# Your paths
sys.path.append("/qfs/projects/gridpack_wind/grid2op_interface/grid2op_local/")
sys.path.append("/qfs/projects/gridpack_wind/grid2op_interface/GridPACK/python/src/")

import grid2op
from grid2op import make
from grid2op.Action import BaseAction
from grid2op.Reward import BaseReward
from grid2op.gym_compat import GymEnv, BoxGymObsSpace, DiscreteActSpace

import stable_baselines3
from stable_baselines3 import DQN
from stable_baselines3.common.vec_env import DummyVecEnv, VecFrameStack
from stable_baselines3.common.callbacks import BaseCallback

import gymnasium as gym
from grid2op_backend import GridPACKBackend

# -------------------- Constants & Defaults --------------------
T_FAULT_CLEAR = 2.0
T_FINAL_CHECK = T_FAULT_CLEAR + 4.0
LOAD_SHED_RATIO = 0.2
C1, C2, C3 = 1.0, 10.0, 100.0

DEFAULT_GRID_XML = "input_9bus.xml"
DEFAULT_CONFIG_PATH = "/qfs/projects/gridpack_wind/grid2op_interface/GridPACK/python/src/grid2op_integration/test_grid2op"
DEFAULT_CONTROLLED_LOADS = [4, 5]
OBS_ATTRS = ["load_v", "load_p"]  # extend if you want

# -------------------- Reward --------------------
class FIDVRReward(BaseReward):
    def __init__(self):
        super().__init__()
        self.init_loads = None
        self.previous_loads = None

    def __call__(self, action, env, has_error, is_done, is_illegal, is_ambiguous) -> float:
        obs = env.get_obs()
        if obs is None: return 0.0

        if self.init_loads is None:
            self.init_loads = obs.load_p.copy()
            self.previous_loads = self.init_loads.copy()

        v = obs.load_v
        p = obs.load_p
        t = obs.current_step * env._time_step

        if t == 0.0:
            self.previous_loads = self.init_loads.copy()
            return 0.0

        if t > T_FINAL_CHECK and (v < 0.95).any():
            return -1000.0

        v_pen = 0.0
        if t > T_FAULT_CLEAR:
            dt_ = t - T_FAULT_CLEAR
            thr = 0.7 if dt_ < 0.33 else 0.8 if dt_ < 0.5 else 0.9 if dt_ < 1.5 else 0.95
            v_pen += np.sum(np.minimum(v - thr, 0.0))

        eps = 1e-6
        load_diff = (self.previous_loads - p) / np.maximum(self.init_loads, eps)
        shed_pen = np.sum(np.maximum(load_diff, 0.0))
        self.previous_loads = p.copy()

        invalid = (p < 1e-3) & (load_diff > 1e-6)
        invalid_pen = float(np.sum(invalid))

        return float(C1 * v_pen - C2 * shed_pen - C3 * invalid_pen)

# -------------------- Standalone Logger (no Gym.Wrapper) --------------------
class TrajectoryLogger:
    """
    Pure-Python logger: call begin_episode() once, then log_step() each step,
    then end_episode() when done. Writes CSV + JSON.
    """
    def __init__(self, outdir: str, attrs: List[str], time_step: float, n_stack: int = 1):
        self.outdir = outdir
        self.attrs = attrs
        self.dt = float(time_step)
        self.n_stack = n_stack
        self.episodes_dir = os.path.join(self.outdir, "episodes")
        os.makedirs(self.episodes_dir, exist_ok=True)
        self.rows = []
        self.ep_id = None
        self.sizes = None
        self.t0 = None

    def _compute_sizes(self, env_native_obs):
        return [getattr(env_native_obs, a).size for a in self.attrs]

    def _split_last_frame(self, obs_vec):
        """
        obs_vec can be stacked (VecFrameStack). We log the last frame only.
        """
        vec = np.asarray(obs_vec).ravel()
        if self.n_stack > 1:
            per = vec.size // self.n_stack
            vec = vec[-per:]  # last frame
        parts, start = [], 0
        for k in self.sizes:
            parts.append(vec[start:start+k])
            start += k
        flat = {}
        for a, part in zip(self.attrs, parts):
            for i, v in enumerate(part):
                flat[f"{a}[{i}]"] = float(v)
        return flat

    def begin_episode(self, initial_obs_vec, env_native_obs):
        self.rows = []
        self.ep_id = uuid.uuid4().hex[:8]
        self.t0 = time.time()
        self.sizes = self._compute_sizes(env_native_obs)
        # initial row (action=-1)
        flat = self._split_last_frame(initial_obs_vec)
        row = dict(step=0, t=0.0, action_idx=-1, reward=0.0, done=False)
        row.update(flat)
        self.rows.append(row)

    def log_step(self, step_idx: int, action_idx: int, reward: float, done: bool, new_obs_vec):
        flat = self._split_last_frame(new_obs_vec)
        row = dict(step=step_idx, t=step_idx * self.dt,
                   action_idx=int(action_idx), reward=float(reward), done=bool(done))
        row.update(flat)
        self.rows.append(row)

    def end_episode(self, info: Dict[str, Any]):
        csv_path = os.path.join(self.episodes_dir, f"ep_{self.ep_id}.csv")
        info_path = os.path.join(self.episodes_dir, f"ep_{self.ep_id}_info.json")
        # CSV
        if self.rows:
            keys = list(self.rows[0].keys())
            with open(csv_path, "w") as f:
                f.write(",".join(keys) + "\n")
                for r in self.rows:
                    f.write(",".join(str(r[k]) for k in keys) + "\n")
        # JSON
        ep_return = float(sum(r["reward"] for r in self.rows))
        ep_len = int(self.rows[-1]["step"])
        meta = dict(
            episode_id=self.ep_id,
            episode_length=ep_len,
            episode_return=ep_return,
            wall_seconds=time.time() - self.t0,
            terminal_info=info or {}
        )
        with open(info_path, "w") as f:
            json.dump(meta, f, indent=2)
        # reset
        self.rows, self.ep_id, self.sizes, self.t0 = [], None, None, None

# -------------------- Env + Actions --------------------
def make_env(config_path: str, grid_xml: str, time_step_s: float) -> GymEnv:
    # read XML to get timestep
    with open(grid_xml, 'r') as f:
        data = f.read()
    bs_data = BeautifulSoup(data, features="lxml")
    gridpack_stepsize = float(bs_data.find('timestep').text)
    
    print(config_path, grid_xml)
    env = make(
        config_path,
        grid_path=grid_xml,
        backend=GridPACKBackend(
            log_freq=1,
            gridpack_stepsize=gridpack_stepsize,
            grid2op_stepsize=time_step_s,  # internal integration step
            can_be_copied=True,
            grid_path=grid_xml
        ),
        data_feeding_kwargs={"time_interval": timedelta(seconds=time_step_s)},
        reward_class=FIDVRReward,
        action_class=BaseAction
    )
    # breakpoint()
    gym_env = GymEnv(env)
    gym_env.observation_space = BoxGymObsSpace(
        env.observation_space,
        attr_to_keep=OBS_ATTRS
    )
    return env, gym_env

def build_action_list(env, controlled_loads: List[int], shed_ratio: float):
    init = env.get_obs()
    p0, q0 = init.load_p.copy(), init.load_q.copy()
    acts, n = [], len(controlled_loads)
    for mask in range(2 ** n):
        new_p, new_q = p0.copy(), q0.copy()
        for bit, lid in enumerate(controlled_loads):
            if (mask >> bit) & 1:
                frac = 1.0 - float(shed_ratio)
                new_p[lid] = frac * p0[lid]
                new_q[lid] = frac * q0[lid]
        acts.append(env.action_space({"injection": {"load_p": new_p, "load_q": new_q}}))
    return acts

# -------------------- Config --------------------
@dataclass
class RunConfig:
    mode: str                   # "train" | "infer-none" | "infer-agent"
    sim_seconds: float
    time_step_s: float
    n_stack: int
    episodes: int
    seed: Optional[int]
    config_path: str
    grid_xml: str
    outdir: str
    model_path: str
    controlled_loads: List[int]
    shed_ratio: float
    # training
    total_timesteps: int = 20_000
    learning_rate: float = 1e-3
    buffer_size: int = 10000
    batch_size: int = 64
    learning_starts: int = 1000
    train_freq: int = 4
    target_update_interval: int = 500
    verbose: int = 1

def save_run_meta(cfg: RunConfig):
    meta = {
        "run_id": uuid.uuid4().hex[:8],
        "timestamp": dt.datetime.utcnow().isoformat()+"Z",
        "grid2op_version": grid2op.__version__,
        "stable_baselines3_version": stable_baselines3.__version__,
        "gymnasium_version": gym.__version__,
        "config": asdict(cfg)
    }
    os.makedirs(cfg.outdir, exist_ok=True)
    with open(os.path.join(cfg.outdir, "meta.json"), "w") as f:
        json.dump(meta, f, indent=2)

# -------------------- Training logger (SB3 callback, no Gym.Wrapper) --------------------
class TrainingLoggerCallback(BaseCallback):
    """
    Logs transitions from the SB3 training loop.
    Assumes single env (n_envs=1). No Gym.Wrapper used.
    """
    def __init__(self, outdir: str, attrs: List[str], time_step: float, n_stack: int):
        super().__init__()
        self.logger = TrajectoryLogger(outdir, attrs, time_step, n_stack)
        self.step_idx = 0
        self.episode_active = False

    def _on_rollout_start(self) -> None:
        # Start of a rollout: begin episode if not active (use current obs)
        try:
            # current observation before steps:
            obs = self.model._last_obs  # (n_envs, obs_dim)
            env_native = self.model.env.envs[0].env.get_obs()
            self.logger.begin_episode(obs[0], env_native)
            self.step_idx = 0
            self.episode_active = True
        except Exception:
            pass

    def _on_step(self) -> bool:
        # SB3 exposes these in self.locals for off-policy algos
        obs = self.locals.get("new_obs", None)   # (n_envs, obs_dim)
        actions = self.locals.get("actions", None)  # (n_envs,)
        rewards = self.locals.get("rewards", None)  # (n_envs,)
        dones = self.locals.get("dones", None)      # (n_envs,)
        infos = self.locals.get("infos", None)      # list of dicts

        if obs is None or actions is None:  # nothing to log
            return True

        a = int(np.asarray(actions)[0])
        r = float(np.asarray(rewards)[0]) if rewards is not None else 0.0
        d = bool(np.asarray(dones)[0]) if dones is not None else False
        self.step_idx += 1
        self.logger.log_step(self.step_idx, a, r, d, obs[0])

        if d:
            info0 = infos[0] if isinstance(infos, (list, tuple)) and infos else {}
            self.logger.end_episode(info0)
            self.episode_active = False

        return True

    def _on_training_end(self) -> None:
        # If training stops mid-episode, flush whatever we have
        if self.episode_active and self.logger.rows:
            self.logger.end_episode({"note": "training_end_flush"})

# -------------------- Inference (manual loops) --------------------
def run_infer_none(cfg: RunConfig):
    env, gym_env = make_env(cfg.config_path, cfg.grid_xml, cfg.time_step_s)
    gym_env.action_space = DiscreteActSpace(env.action_space,
                                         action_list=build_action_list(env, cfg.controlled_loads, cfg.shed_ratio))
    logger = TrajectoryLogger(cfg.outdir, OBS_ATTRS, cfg.time_step_s, n_stack=1)

    for ep in range(cfg.episodes):
        print(f"[infer-none] episode {ep+1}/{cfg.episodes}")
        obs, info = gym_env.reset()
        logger.begin_episode(obs, env.get_obs())
        max_steps = int(math.ceil(cfg.sim_seconds / cfg.time_step_s))
        done = False
        for k in range(1, max_steps+1):
            breakpoint()
            obs, reward, terminated, truncated, info = gym_env.step(0)  # noop = index 0
            done = bool(terminated or truncated)
            logger.log_step(k, 0, reward, done, obs)
            if done: break
        logger.end_episode(info or {})
    gym_env.close()

def run_infer_agent(cfg: RunConfig):
    env, gym_env = make_env(cfg.config_path, cfg.grid_xml, cfg.time_step_s)
    gym_env.action_space = DiscreteActSpace(env.action_space,
                                         action_list=build_action_list(env, cfg.controlled_loads, cfg.shed_ratio))
    venv = DummyVecEnv([lambda: gym_env])
    venv = VecFrameStack(venv, n_stack=cfg.n_stack)
    model = DQN.load(cfg.model_path)

    logger = TrajectoryLogger(cfg.outdir, OBS_ATTRS, cfg.time_step_s, n_stack=cfg.n_stack)

    for ep in range(cfg.episodes):
        print(f"[infer-agent] episode {ep+1}/{cfg.episodes}")
        obs = venv.reset()
        logger.begin_episode(obs[0], env.get_obs())
        max_steps = int(math.ceil(cfg.sim_seconds / cfg.time_step_s))
        for k in range(1, max_steps+1):
            action, _ = model.predict(obs, deterministic=True)
            obs, reward, term, trunc, info = venv.step(action)
            done = bool(term[0] or trunc[0])
            logger.log_step(k, int(action[0]), float(reward[0]), done, obs[0])
            if done: break
        logger.end_episode(info[0] if isinstance(info, (list, tuple)) and info else {})
    venv.close()

# -------------------- Training --------------------
def run_train(cfg: RunConfig):
    env, gym_env = make_env(cfg.config_path, cfg.grid_xml, cfg.time_step_s)
    gym_env.action_space = DiscreteActSpace(env.action_space,
                                         action_list=build_action_list(env, cfg.controlled_loads, cfg.shed_ratio))
    venv = DummyVecEnv([lambda: gym_env])
    venv = VecFrameStack(venv, n_stack=cfg.n_stack)

    model = DQN(
        "MlpPolicy",
        venv,
        learning_rate=cfg.learning_rate,
        buffer_size=cfg.buffer_size,
        batch_size=cfg.batch_size,
        learning_starts=cfg.learning_starts,
        train_freq=cfg.train_freq,
        target_update_interval=cfg.target_update_interval,
        verbose=cfg.verbose,
        tensorboard_log=os.path.join(cfg.outdir, "tb")
    )

    cb = TrainingLoggerCallback(cfg.outdir, OBS_ATTRS, cfg.time_step_s, cfg.n_stack)
    model.learn(total_timesteps=cfg.total_timesteps, callback=cb, log_interval=10)
    os.makedirs(os.path.dirname(cfg.model_path), exist_ok=True)
    model.save(cfg.model_path)
    print(f"[train] model saved -> {cfg.model_path}")

# -------------------- CLI --------------------
def parse_args() -> RunConfig:
    p = argparse.ArgumentParser()
    p.add_argument("--mode", choices=["train","infer-none","infer-agent"], default="infer-none")
    p.add_argument("--sim-seconds", type=float, default=300.0)
    p.add_argument("--time-step", type=float, default=0.05)
    p.add_argument("--n-stack", type=int, default=4)
    p.add_argument("--episodes", type=int, default=1)
    p.add_argument("--seed", type=int, default=None)
    p.add_argument("--config-path", type=str, default=DEFAULT_CONFIG_PATH)
    p.add_argument("--grid-xml", type=str, default=DEFAULT_GRID_XML)
    p.add_argument("--outdir", type=str, default=".")
    p.add_argument("--model-path", type=str, default="./models/dqn_fidvr_load_shedding")
    p.add_argument("--controlled-loads", type=str, default="0,1,2")
    p.add_argument("--shed-ratio", type=float, default=LOAD_SHED_RATIO)
    # training
    p.add_argument("--total-timesteps", type=int, default=20000)
    p.add_argument("--learning-rate", type=float, default=1e-3)
    p.add_argument("--buffer-size", type=int, default=10000)
    p.add_argument("--batch-size", type=int, default=64)
    p.add_argument("--learning-starts", type=int, default=1000)
    p.add_argument("--train-freq", type=int, default=4)
    p.add_argument("--target-update-interval", type=int, default=500)
    a = p.parse_args()
    controlled = [int(x) for x in a.controlled_loads.split(",") if x.strip()!=""]
    return RunConfig(
        mode=a.mode,
        sim_seconds=a.sim_seconds,
        time_step_s=a.time_step,
        n_stack=a.n_stack,
        episodes=a.episodes,
        seed=a.seed,
        config_path=a.config_path,
        grid_xml=a.grid_xml,
        outdir=a.outdir,
        model_path=a.model_path,
        controlled_loads=controlled,
        shed_ratio=a.shed_ratio,
        total_timesteps=a.total_timesteps,
        learning_rate=a.learning_rate,
        buffer_size=a.buffer_size,
        batch_size=a.batch_size,
        learning_starts=a.learning_starts,
        train_freq=a.train_freq,
        target_update_interval=a.target_update_interval
    )

def main():
    print(f"Grid2Op: {grid2op.__version__} | SB3: {stable_baselines3.__version__} | Gymnasium: {gym.__version__}")
    cfg = parse_args()
    save_run_meta(cfg)

    if cfg.mode == "infer-none":
        run_infer_none(cfg)
    elif cfg.mode == "infer-agent":
        if not os.path.exists(cfg.model_path + ".zip"):
            raise FileNotFoundError(f"Missing model: {cfg.model_path}.zip")
        run_infer_agent(cfg)
    elif cfg.mode == "train":
        run_train(cfg)
    else:
        raise ValueError(cfg.mode)

if __name__ == "__main__":
    main()
