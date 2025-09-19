import numpy as np
from typing import Dict

import grid2op
from grid2op.Reward import BaseReward
from grid2op.Action import BaseAction
from grid2op.Environment import BaseEnv

from config import T_FINAL_CHECK, T_FAULT_CLEAR, C1, C2, C3

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
        thr = None
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
    
def reward_components(obs, prev_loads, init_loads) -> Dict[str, float]:
    """Mirror FIDVRReward math so we can log parts each step."""
    t = getattr(obs, "current_step", 0) * getattr(obs, "_time_step", 1)
    v = np.asarray(getattr(obs, "load_v")).ravel()
    p = np.asarray(getattr(obs, "load_p")).ravel()
    result = {
        "threshold": np.nan, "v_pen": 0.0, "shed_pen": 0.0, "invalid_pen": 0.0,
        "t_since_fault": 0.0, "fault_cleared": False,
    }
    if t <= 0.0:
        return result
    if t > T_FINAL_CHECK and (v < 0.95).any():
        # We still log the pieces as-is; the final -1000 is applied by reward class.
        pass
    v_pen = 0.0
    thr = None
    if t > T_FAULT_CLEAR:
        dt_ = t - T_FAULT_CLEAR
        thr = 0.7 if dt_ < 0.33 else 0.8 if dt_ < 0.5 else 0.9 if dt_ < 1.5 else 0.95
        v_pen = float(np.sum(np.minimum(v - thr, 0.0)))
        result["threshold"] = float(thr)
        result["t_since_fault"] = float(dt_)
        result["fault_cleared"] = True
    eps = 1e-6
    load_diff = (prev_loads - p) / np.maximum(init_loads, eps)
    shed_pen = float(np.sum(np.maximum(load_diff, 0.0)))
    invalid = (p < 1e-3) & (load_diff > 1e-6)
    invalid_pen = float(np.sum(invalid))

    result["v_pen"] = float(v_pen)
    result["shed_pen"] = float(shed_pen)
    result["invalid_pen"] = float(invalid_pen)
    return result
