from typing import List
import random, collections

import torch
import torch.nn as nn
import torch.optim as optim

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

def build_action_list(env, controlled_loads: List[int], shed_ratio: float):
    init = env.reset()
    p0, q0 = init.load_p.copy(), init.load_q.copy()
    acts, n = [], len(controlled_loads)
    # print(p0, q0)
    # print(env.backend._grid.load)
    # print(env.load_info)
    for mask in range(2 ** n):
        new_p, new_q = p0.copy(), q0.copy()
        for bit, lid in enumerate(controlled_loads):
            # Map "bus id" -> load index used by backend
            # breakpoint()
            lid_processed = env.backend._grid.load.loc[env.backend._grid.load["bus"]==lid].index.values[0]
            if (mask >> bit) & 1:
                frac = 1.0 - float(shed_ratio)
                new_p[lid_processed] = frac * p0[lid_processed]
                new_q[lid_processed] = frac * q0[lid_processed]
        acts.append(env.action_space({"injection": {"load_p": new_p, "load_q": new_q}}))
    return acts