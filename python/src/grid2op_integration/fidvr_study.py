import os, sys, math, random, collections, time
from typing import List, Tuple
import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from datetime import timedelta
from bs4 import BeautifulSoup

sys.path += [
    "/qfs/projects/gridpack_wind/grid2op_interface/grid2op_local/",
    "/qfs/projects/gridpack_wind/grid2op_interface/GridPACK/python/src/",
]

import grid2op
from grid2op import make
from grid2op.Action import BaseAction
from grid2op.Reward import BaseReward
from grid2op.Runner import Runner

from grid2op_backend import GridPACKBackend

# ------------ Config ------------
GRID_XML = "input_9bus.xml"
CONFIG_PATH = "/qfs/projects/gridpack_wind/grid2op_interface/GridPACK/python/src/grid2op_integration/test_grid2op"
OBS_ATTRS = ["load_v", "load_p"]
CONTROLLED_LOADS = [4, 5]
SHED_RATIO = 0.2
TIME_STEP_S = 0.05
TOTAL_STEPS = 20_000
TARGET_SYNC = 1000
BATCH_SIZE = 64
GAMMA = 0.99
LR = 1e-3
EPS_START, EPS_END, EPS_DECAY = 1.0, 0.05, 15_000

# ------------ Reward (as you had) ------------
T_FAULT_CLEAR = 2.0
T_FINAL_CHECK = T_FAULT_CLEAR + 4.0
C1, C2, C3 = 1.0, 10.0, 100.0

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
        v, p = obs.load_v, obs.load_p
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

# ------------ Helpers ------------
def read_gridpack_stepsize(xml_path: str) -> float:
    with open(xml_path, "r") as f: txt = f.read()
    bs = BeautifulSoup(txt, "lxml")
    return float(bs.find("timestep").text)

def build_action_list(env, controlled_loads: List[int], shed_ratio: float):
    init = env.reset()
    p0, q0 = init.load_p.copy(), init.load_q.copy()
    acts, n = [], len(controlled_loads)
    
    # breakpoint()
    for mask in range(2 ** n):
        new_p, new_q = p0.copy(), q0.copy()
        for bit, lid in enumerate(controlled_loads):
            lid_processed = env.backend._grid.load.loc[env.backend._grid.load["bus"]==lid].index.values[0]
            print(bit, lid, lid_processed, p0, q0)
            if (mask >> bit) & 1:
                frac = 1.0 - float(shed_ratio)
                new_p[lid_processed] = frac * p0[lid_processed]
                new_q[lid_processed] = frac * q0[lid_processed]
        acts.append(env.action_space({"injection": {"load_p": new_p, "load_q": new_q}}))
    return acts

def flatten_obs(env_obs) -> np.ndarray:
    parts = [np.asarray(getattr(env_obs, a)).ravel() for a in OBS_ATTRS]
    return np.concatenate(parts).astype(np.float32, copy=False)

# ------------ Q-network ------------
class MLP(nn.Module):
    def __init__(self, in_dim: int, out_dim: int):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(in_dim, 256), nn.ReLU(),
            nn.Linear(256, 256), nn.ReLU(),
            nn.Linear(256, out_dim)
        )
    def forward(self, x): return self.net(x)

# ------------ Replay buffer ------------
Transition = collections.namedtuple("Transition", "s a r s2 done")
class Replay:
    def __init__(self, cap=100000): self.buf, self.cap = [], cap
    def push(self, *args):
        if len(self.buf) >= self.cap: self.buf.pop(0)
        self.buf.append(Transition(*args))
    def sample(self, n):
        batch = random.sample(self.buf, n)
        return Transition(*zip(*batch))
    def __len__(self): return len(self.buf)

# ------------ Build env (pure grid2op) ------------
def make_env():
    stepsize_gp = read_gridpack_stepsize(GRID_XML)
    env = make(
        CONFIG_PATH,
        grid_path=GRID_XML,
        backend=GridPACKBackend(
            log_freq=1,
            gridpack_stepsize=stepsize_gp,
            grid2op_stepsize=TIME_STEP_S,
            can_be_copied=True,
            grid_path=GRID_XML
        ),
        data_feeding_kwargs={"time_interval": timedelta(seconds=TIME_STEP_S)},
        reward_class=FIDVRReward,
        action_class=BaseAction
    )
    return env

# ------------ Training loop (no Gym) ------------
def train():
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    env = make_env()
    action_list = build_action_list(env, CONTROLLED_LOADS, SHED_RATIO)

    # obs dims
    o0, _ = env.reset()
    s0 = flatten_obs(o0)
    in_dim, out_dim = s0.size, len(action_list)

    q = MLP(in_dim, out_dim).to(device)
    qt = MLP(in_dim, out_dim).to(device)
    qt.load_state_dict(q.state_dict())
    opt = optim.Adam(q.parameters(), lr=LR)
    replay = Replay(100000)

    global_step = 0
    ep_idx = 0
    best_return = -1e9

    def epsilon(step):
        # linear decay
        if step >= EPS_DECAY: return EPS_END
        return EPS_START + (EPS_END - EPS_START) * (step / EPS_DECAY)

    while global_step < TOTAL_STEPS:
        ep_idx += 1
        obs, info = env.reset()
        s = flatten_obs(obs)
        ep_ret = 0.0
        done = False

        # optional cap: simulate ~300s
        max_steps = int(math.ceil(300.0 / TIME_STEP_S))

        for t in range(max_steps):
            eps = epsilon(global_step)
            if random.random() < eps:
                a_idx = random.randrange(out_dim)
            else:
                with torch.no_grad():
                    qs = q(torch.from_numpy(s).unsqueeze(0).to(device))
                    a_idx = int(torch.argmax(qs, dim=1).item())

            # step env natively (no Gym)
            next_obs, r, term, trunc, info = env.step(action_list[a_idx])
            d = bool(term or trunc)
            s2 = flatten_obs(next_obs)
            replay.push(s, a_idx, float(r), s2, d)
            s = s2
            ep_ret += float(r)
            global_step += 1

            # learn
            if len(replay) >= BATCH_SIZE:
                batch = replay.sample(BATCH_SIZE)
                ss = torch.tensor(np.stack(batch.s), dtype=torch.float32, device=device)
                aa = torch.tensor(batch.a, dtype=torch.int64, device=device)
                rr = torch.tensor(batch.r, dtype=torch.float32, device=device)
                ss2 = torch.tensor(np.stack(batch.s2), dtype=torch.float32, device=device)
                dd = torch.tensor(batch.done, dtype=torch.float32, device=device)

                qsa = q(ss).gather(1, aa.unsqueeze(1)).squeeze(1)
                with torch.no_grad():
                    maxq = qt(ss2).max(1).values
                    target = rr + (1.0 - dd) * GAMMA * maxq
                loss = nn.functional.smooth_l1_loss(qsa, target)

                opt.zero_grad()
                loss.backward()
                nn.utils.clip_grad_norm_(q.parameters(), 1.0)
                opt.step()

            # target sync
            if global_step % TARGET_SYNC == 0:
                qt.load_state_dict(q.state_dict())

            if d or global_step >= TOTAL_STEPS:
                print(f"ep {ep_idx:04d} | steps {global_step} | return {ep_ret:.3f} | eps {eps:.3f}")
                if ep_ret > best_return:
                    best_return = ep_ret
                    os.makedirs("./models", exist_ok=True)
                    torch.save(q.state_dict(), "./models/dqn_grid2op.pt")
                break

    env.close()
    print("Training done.")

if __name__ == "__main__":
    print(f"Grid2Op {grid2op.__version__} | Torch {torch.__version__}")
    train()
