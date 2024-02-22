def test_probe_run_gym(
    gym_env,
    CONTROLLED_LOADS,
    max_steps=60,
    T_FAULT_CLEAR=2.0,
    ts=0.05,
    n_stack=4,
    obs_order="vp",  # "vp" if per-frame is [load_v, load_p]; use "pv" if it's [load_p, load_v]
):
    """
    Probe the SAME gym_env used for training (GymEnv or DummyVecEnv+VecFrameStack).
    Prints time, TVRC phase, min(load_v), controlled loads (MW), per-unit shed, reward,
    and whether (t > T_pf + 4s AND minV < 0.95) holds.

    Robust to:
      - (obs, info) tuples on reset/step
      - 4-tuple vs 5-tuple step APIs
      - Stacked obs being flat (feat*n_stack,) OR 2-D (n_stack, feat)
    """
    import numpy as np

    T_FINAL_CHECK = T_FAULT_CLEAR + 4.0
    is_vec = hasattr(gym_env, "num_envs")

    def tvrc_phase(t):
        if t < T_FAULT_CLEAR:
            return "pre-clear"
        dt = t - T_FAULT_CLEAR
        if dt < 0.33:  return "TVRC[0.7] (0–0.33s)"
        if dt < 0.5:   return "TVRC[0.8] (0.33–0.5s)"
        if dt < 1.5:   return "TVRC[0.9] (0.5–1.5s)"
        return "TVRC[0.95] (>1.5s)"

    def _reset_obs(env):
        out = env.reset()
        obs = out[0] if isinstance(out, tuple) and len(out) == 2 else out
        if is_vec and hasattr(obs, "ndim") and obs.ndim >= 2:
            obs = obs[0]  # env 0
        return obs  # shape: (feat*n_stack,) or (n_stack, feat) or (feat,)

    def _step_any(env, a_idx):
        if is_vec:
            obs, rewards, dones, infos = env.step([a_idx])
            r = float(np.array(rewards).ravel()[0])
            d = bool(np.array(dones).ravel()[0])
            o = obs[0] if hasattr(obs, "ndim") and obs.ndim >= 2 else obs
            return o, r, d, infos
        else:
            out = env.step(a_idx)
            if isinstance(out, tuple):
                if len(out) == 5:
                    o, r, term, trunc, inf = out
                    d = bool(term or trunc)
                elif len(out) == 4:
                    o, r, d, inf = out
                else:
                    # best-effort fallback
                    o, r, d, inf = out[0], 0.0, False, {}
            else:
                o, r, d, inf = out, 0.0, False, {}
            return o, float(r), bool(d), inf

    def _last_frame(obs_vec):
        """Return 1-D last frame vector; handles (feat*n_stack,), (n_stack, feat), (feat,)"""
        arr = np.asarray(obs_vec)
        if arr.ndim == 1:
            # flat; try to split by n_stack
            if n_stack > 1 and arr.size % n_stack == 0:
                per = arr.size // n_stack
                return arr[-per:]
            return arr  # unstacked
        elif arr.ndim == 2:
            # stacked first
            if arr.shape[0] == n_stack:
                return arr[-1]
            # unexpected 2-D, flatten last slice
            return arr.reshape(-1)[-arr.shape[-1]:]
        else:
            # higher dims: flatten
            return arr.reshape(-1)

    def _split_frame(frame_vec):
        """Split frame into (load_v, load_p) per obs_order."""
        v = np.asarray(frame_vec)
        assert v.ndim == 1, "Frame must be 1-D after _last_frame."
        if v.size % 2 != 0:
            raise ValueError(f"Frame length {v.size} not even; expected [load_v, load_p] or [load_p, load_v].")
        half = v.size // 2
        if obs_order.lower() == "vp":
            load_v = v[:half]
            load_p = v[half:]
        else:  # "pv"
            load_p = v[:half]
            load_v = v[half:]
        return load_v.astype(float), load_p.astype(float)

    # ==== init ====
    obs0 = _reset_obs(gym_env)
    frame0 = _last_frame(obs0)
    load_v0, load_p0 = _split_frame(frame0)
    init_loads = load_p0.copy()
    eps = 1e-6

    print("\n=== Probe run (first ~{} steps) ===".format(max_steps))
    print("step |   t[s] | phase                | minV(loads) |  loads(controlled) [MW]  | shed_pu(controlled) |   reward   | -1000?")
    print("-"*125)

    NOOP_IDX = 0
    ALL_ONES_IDX = (1 << len(CONTROLLED_LOADS)) - 1
    already_shed = False
    cum_reward = 0.0
    done_flag = False

    for step in range(max_steps):
        t = step * ts

        a_idx = ALL_ONES_IDX if (t >= T_FAULT_CLEAR + ts and not already_shed) else NOOP_IDX
        already_shed = already_shed or (a_idx == ALL_ONES_IDX)

        obs_out, r, d, info = _step_any(gym_env, a_idx)
        frame = _last_frame(obs_out)
        load_v, load_p = _split_frame(frame)

        min_v = float(np.min(load_v))
        shed_pu = np.maximum((init_loads - load_p) / np.maximum(init_loads, eps), 0.0)

        # guard against index range mismatches
        max_idx = len(load_p) - 1
        bad = [i for i in CONTROLLED_LOADS if i > max_idx]
        if bad:
            ctrl_idxs = [i for i in CONTROLLED_LOADS if i <= max_idx]
        else:
            ctrl_idxs = CONTROLLED_LOADS

        loads_now_ctrl = " ".join([f"{load_p[i]:6.1f}" for i in ctrl_idxs]) if ctrl_idxs else "(no valid idx)"
        shed_ctrl      = " ".join([f"{shed_pu[i]:4.2f}" for i in ctrl_idxs]) if ctrl_idxs else "(no valid idx)"

        term_flag = (t > T_FINAL_CHECK) and (min_v < 0.95)
        cum_reward += r

        if step % 5 == 0 or d:
            print(f"{step:4d} | {t:6.2f} | {tvrc_phase(t):20s} |     {min_v:6.3f} | [{loads_now_ctrl}] | [{shed_ctrl}] | {r:9.3f} | {'YES' if term_flag else 'no '}")

        if d:
            done_flag = True
            print(f"\nEpisode ended at step {step} (t={t:.2f}s). Cumulative return: {cum_reward:.3f}")
            break

    if not done_flag:
        print(f"\nStopped after {max_steps} steps. Cumulative return: {cum_reward:.3f}")
