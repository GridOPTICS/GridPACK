import os, sys, time, math, psutil, random

sys.path += [
    "/qfs/projects/gridpack_wind/grid2op_interface/grid2op_local/",
    "/qfs/projects/gridpack_wind/grid2op_interface/GridPACK/python/src/",
]

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.optim as optim
from datetime import timedelta

import pathlib
import argparse

from config import CONFIG_DIR, TIME_STEP_S, LR, GAMMA, TARGET_SYNC, BATCH_SIZE, EPS_START, EPS_END, EPS_DECAY, CONTROLLED_LOADS, SHED_RATIO, TOTAL_STEPS, SAVE_EVERY_STEPS

from model import MLP, Replay, build_action_list

from data_loggers import TrainLogger, SimDataSaver, TimingLogger
from reward import FIDVRReward, reward_components
from utils import read_gridpack_stepsize, flatten_obs, get_sizes, detect_attrs



import grid2op
from grid2op import make
from grid2op.Action import BaseAction
from grid2op.Reward import BaseReward
from grid2op.Runner import Runner

from grid2op_backend import GridPACKBackend

# ------------ Build env (pure grid2op) ------------
def make_env(grid_xml, config_path):
    gridpack_info = read_gridpack_stepsize(grid_xml)
    stepsize_gp = gridpack_info["global_timeStep"]
    fault_clearance = gridpack_info["faultEvents"][0]["endFault"]
    env = make(
        config_path,
        grid_path=grid_xml,
        backend=GridPACKBackend(
            log_freq=1,
            gridpack_stepsize=stepsize_gp,
            grid2op_stepsize=TIME_STEP_S,
            can_be_copied=True,
            grid_path=grid_xml
        ),
        data_feeding_kwargs={"time_interval": timedelta(seconds=TIME_STEP_S)},
        reward_class=FIDVRReward(fault_clearance),
        action_class=BaseAction
    )
    return env, fault_clearance

def run_simulation(grid_xml: str, config_path: str, outdir: str, episodes: int, max_steps: int, log_each_step: bool=True):
    """Do-nothing agent (no shedding), just run and log."""
    timelog = TimingLogger(outdir=outdir, run_name="simulate", enable_traces=True)
    
    timelog.start("setup_s")
    env, t_fault_clear = make_env(grid_xml, config_path)
    # Prepare logging based on first obs
    obs0 = env.reset()
    
    attrs_present = detect_attrs(obs0)
    n_bus, n_load = get_sizes(obs0)
    
    # Build an action list but use index 0 as "do nothing"
    actions = build_action_list(env, CONTROLLED_LOADS[grid_xml], SHED_RATIO)
    do_nothing_idx = 0
    timelog.stop("setup_s")
    
    timelog.start("logging_s")
    saver = SimDataSaver(
        env.backend, 
        outdir, 
        run_name="simulate", 
        n_bus=n_bus, 
        n_load=n_load, 
        attrs_present=attrs_present
    )
    timelog.stop("logging_s")

    timelog.start("simulate_s")
    proc = psutil.Process(os.getpid())

    for ep in range(1, episodes+1):
        print("============= Starting Episode:", ep, "==============")
        timelog.begin_episode(ep)

        obs = env.reset()
        
        # init for reward components
        init_loads = np.asarray(getattr(obs, "load_p")).copy()
        prev_loads = init_loads.copy()

        for t in range(max_steps):
            a_idx = do_nothing_idx
            
            t0 = time.perf_counter()
            next_obs, r, done, info = env.step(actions[a_idx])
            step_dt_ms = (time.perf_counter() - t0) * 1000.0
            rss_mib = proc.memory_info().rss / (1024.0 * 1024.0)
            timelog.add_step(step_dt_ms, rss_mib)
            
            # log arrays
            
            if log_each_step:
                tlog0 = time.perf_counter()
                step_time = env.backend._counter_time
                rc = reward_components(env, step_time, t_fault_clear, next_obs, prev_loads, init_loads)
                saver.log_arrays(next_obs)
                saver.log_step(
                    env.backend, 
                    episode=ep, 
                    t=next_obs.current_step, 
                    step=t, 
                    a_idx=a_idx, 
                    reward=float(r), 
                    rc=rc
                )
                timelog.sections["logging_s"] += (time.perf_counter() - tlog0)

            prev_loads = np.asarray(getattr(next_obs, "load_p")).copy()
            if done: break

        ep_meta = {
            "system": pathlib.Path(grid_xml).stem,  # e.g., "case9" / "IEEE39bus"
            "agent": "DoNothing",
            "logging": "on" if log_each_step else "off",
            "episode_seconds": getattr(next_obs, "current_step", 0)  # fill if you prefer
        }
        timelog.end_episode(meta=ep_meta)

    timelog.stop("simulate_s")

    timelog.start("logging_s")
    saver.close()
    env.close()
    timelog.stop("logging_s")
    
    timelog.dump_summaries()
    print(f"[simulate] done. Logs in {pathlib.Path(outdir) / 'simulate'}")

def train(grid_xml: str, config_path: str, outdir: str, max_steps: int, total_steps: int, log_each_step: bool=False):
    timelog = TimingLogger(outdir=outdir, run_name="train", enable_traces=True)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    proc = psutil.Process(os.getpid())
    
    timelog.start("setup_s")
    env, t_fault_clear = make_env(grid_xml, config_path)
    action_list = build_action_list(env, CONTROLLED_LOADS[grid_xml], SHED_RATIO)
    
    o0 = env.reset()
    s0 = flatten_obs(o0)
    in_dim, out_dim = s0.size, len(action_list)

    q = MLP(in_dim, out_dim).to(device)
    qt = MLP(in_dim, out_dim).to(device)
    qt.load_state_dict(q.state_dict())
    opt = optim.Adam(q.parameters(), lr=LR)
    replay = Replay(100000)

    logger = TrainLogger(outdir=str(pathlib.Path(outdir) / "train"))
    
    # Optional sim saver during training
    saver = None
    if log_each_step:
        attrs_present = detect_attrs(o0)
        n_bus, n_load = get_sizes(o0)
        saver = SimDataSaver(
            env.backend,
            outdir=outdir, 
            run_name="train_sim", 
            n_bus=n_bus, 
            n_load=n_load, 
            attrs_present=attrs_present
        )
    timelog.stop("setup_s")

    timelog.start("training_s")
    global_step = 0
    ep_idx = 0
    best_return = -1e9

    def epsilon(step):
        # linear decay
        if step >= EPS_DECAY: return EPS_END
        return EPS_START + (EPS_END - EPS_START) * (step / EPS_DECAY)

    loop_t0 = time.perf_counter()
    while global_step < total_steps:
        print("============= Starting Episode:", ep_idx+1, "==============")
        ep_idx += 1
        timelog.begin_episode(ep_idx)

        obs = env.reset()
        s = flatten_obs(obs)
        ep_ret = 0.0
        done = False
        ep_len = 0
        last_eps = epsilon(global_step)

        # reward component tracking init
        init_loads = np.asarray(getattr(obs, "load_p")).copy()
        prev_loads = init_loads.copy()

        for t in range(max_steps):
            eps = epsilon(global_step)
            last_eps = eps
            if random.random() < eps:
                a_idx = random.randrange(out_dim)
            else:
                with torch.no_grad():
                    qs = q(torch.from_numpy(s).unsqueeze(0).to(device))
                    a_idx = int(torch.argmax(qs, dim=1).item())

            # -------- env.step timing (env latency + RSS)
            t0 = time.perf_counter()
            next_obs, r, d, info = env.step(action_list[a_idx])
            step_dt_ms = (time.perf_counter() - t0) * 1000.0
            rss_mib = proc.memory_info().rss / (1024.0 * 1024.0)
            timelog.add_step(step_dt_ms, rss_mib)

            s2 = flatten_obs(next_obs)
            replay.push(s, a_idx, float(r), s2, d)
            s = s2
            ep_ret += float(r)
            ep_len += 1
            global_step += 1

            # learn
            loss_value = None
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
                loss_value = float(loss.item())

            tlog0 = time.perf_counter()
            logger.log_step(global_step, ep_idx, t, r, a_idx, eps, (loss_value if loss_value is not None else None))

            if saver is not None:
                # breakpoint()
                step_time = env.backend._counter_time
                rc = reward_components(env, step_time, t_fault_clear, next_obs, prev_loads, init_loads)
                saver.log_arrays(next_obs)
                saver.log_step(
                    env.backend, 
                    episode=ep_idx, 
                    t=next_obs.current_step, 
                    step=t, 
                    a_idx=a_idx, 
                    reward=float(r), 
                    rc=rc
                )

            # target sync
            if global_step % TARGET_SYNC == 0:
                qt.load_state_dict(q.state_dict())

            # periodic checkpoint + flush
            if global_step % SAVE_EVERY_STEPS == 0:
                torch.save(q.state_dict(), pathlib.Path(outdir) / "train" / f"ckpts/step_{global_step}.pt")
                logger.flush()
            timelog.sections["logging_s"] += (time.perf_counter() - tlog0)

            prev_loads = np.asarray(getattr(next_obs, "load_p")).copy()

            if d or global_step >= total_steps:
                tlog0 = time.perf_counter()
                print(f"ep {ep_idx:04d} | steps {global_step} | return {ep_ret:.3f} | eps {eps:.3f}")
                if ep_ret > best_return:
                    best_return = ep_ret
                    os.makedirs("./models", exist_ok=True)
                    torch.save(q.state_dict(), "./models/dqn_grid2op.pt")
                logger.log_episode(ep_idx, global_step, ep_ret, ep_len, last_eps)
                logger.flush()

                timelog.sections["logging_s"] += (time.perf_counter() - tlog0)

                # record episode-level system efficiency summary
                ep_meta = {
                    "system": pathlib.Path(grid_xml).stem,
                    "agent": "DQN",
                    "logging": "on" if log_each_step else "off",
                    "episode_seconds": getattr(next_obs, "current_step", 0)
                }
                timelog.end_episode(meta=ep_meta)

                break

    timelog.stop("training_s")

    env.close()
    if saver is not None: saver.close()
    logger.close()
    
    timelog.dump_summaries()
    print("[train] done.")

def _model_num_bytes(module: torch.nn.Module) -> int:
    """Parameter + buffer bytes (dtype-aware)."""
    nbytes = 0
    for t in module.state_dict().values():
        nbytes += t.numel() * t.element_size()
    return nbytes

def deploy(grid_xml: str, config_path: str, model_path: str, outdir: str, episodes: int, max_steps: int):
    timelog = TimingLogger(outdir=outdir, run_name="deploy", enable_traces=True)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    proc = psutil.Process(os.getpid())

    # -------- Setup --------
    timelog.start("setup_s")
    env, t_fault_clear = make_env(grid_xml, config_path)
    actions = build_action_list(env, CONTROLLED_LOADS[grid_xml], SHED_RATIO)

    # Model dims from env
    o0 = env.reset()
    in_dim, out_dim = flatten_obs(o0).size, len(actions)
    q = MLP(in_dim, out_dim).to(device)
    q.load_state_dict(torch.load(model_path, map_location=device))
    q.eval()

    # Static policy footprint
    model_bytes = _model_num_bytes(q)
    model_mib = model_bytes / (1024.0 * 1024.0)

    if device.type == "cuda":
        torch.cuda.reset_peak_memory_stats()
        # force param residency
        _ = sum(p.numel() for p in q.parameters())

    attrs_present = detect_attrs(o0)
    n_bus, n_load = get_sizes(o0)

    saver = SimDataSaver(
        env.backend,
        outdir=outdir,
        run_name="deploy",
        n_bus=n_bus,
        n_load=n_load,
        attrs_present=attrs_present
    )
    timelog.stop("setup_s")

    outdir_path = pathlib.Path(outdir) / "deploy"
    outdir_path.mkdir(parents=True, exist_ok=True)

    # -------- Episodes --------
    for ep in range(1, episodes + 1):
        timelog.begin_episode(ep)

        obs = env.reset()
        s = flatten_obs(obs)
        init_loads = np.asarray(getattr(obs, "load_p")).copy()
        prev_loads = init_loads.copy()

        # Per-episode accumulators for summaries
        infer_ms_acc = []   # forward-pass (ms)
        env_ms_acc = []     # env.step (ms)
        loop_ms_acc = []    # closed-loop (ms)
        rho_vals = []       # real-time factor
        rss_mib_series = []

        # Δt estimate (simulation step)
        sim_dt_hint = getattr(env.backend, "_dt", None)
        prev_sim_t = getattr(obs, "current_step", None)

        for t in range(max_steps):
            loop_t0 = time.perf_counter()

            # ---- Inference
            if device.type == "cuda":
                torch.cuda.synchronize()
            t0 = time.perf_counter()
            with torch.no_grad():
                qs = q(torch.from_numpy(s).unsqueeze(0).to(device))
                a_idx = int(torch.argmax(qs, dim=1).item())
            if device.type == "cuda":
                torch.cuda.synchronize()
            inf_ms = (time.perf_counter() - t0) * 1000.0
            infer_ms_acc.append(inf_ms)
            timelog.sections["inference_s"] += (inf_ms / 1000.0)

            # ---- Environment step
            t0 = time.perf_counter()
            next_obs, r, done, info = env.step(actions[a_idx])
            env_ms = (time.perf_counter() - t0) * 1000.0
            env_ms_acc.append(env_ms)

            # ---- Logging (arrays + reward components)
            tlog0 = time.perf_counter()
            step_time = env.backend._counter_time
            rc = reward_components(env, step_time, t_fault_clear, next_obs, prev_loads, init_loads)
            saver.log_arrays(next_obs)
            saver.log_step(
                env.backend,
                episode=ep,
                t=next_obs.current_step,
                step=t,
                a_idx=a_idx,
                reward=float(r),
                rc=rc
            )
            log_ms = (time.perf_counter() - tlog0) * 1000.0
            timelog.sections["logging_s"] += (log_ms / 1000.0)

            # ---- Closed-loop latency
            loop_ms_val = (time.perf_counter() - loop_t0) * 1000.0
            loop_ms_acc.append(loop_ms_val)

            # ---- Memory + ρ
            rss_mib = proc.memory_info().rss / (1024.0 * 1024.0)
            rss_mib_series.append(rss_mib)

            sim_t = getattr(next_obs, "current_step", None)
            if sim_dt_hint is not None:
                dt_sim = float(sim_dt_hint)
            elif (prev_sim_t is not None) and (sim_t is not None):
                dt_sim = float(sim_t - prev_sim_t)
            else:
                dt_sim = None
            rho = (dt_sim / (loop_ms_val / 1000.0)) if dt_sim is not None else np.nan
            if np.isfinite(rho):
                rho_vals.append(float(rho))
            prev_sim_t = sim_t

            # ---- Emit per-step trace via TimingLogger (extended)
            timelog.add_step(
                latency_ms=env_ms,
                rss_mib=rss_mib,
                infer_ms=inf_ms,
                loop_ms=loop_ms_val,
                log_ms=log_ms,
                rho=(float(rho) if np.isfinite(rho) else np.nan),
                a_idx=int(a_idx),
                t=(float(sim_t) if sim_t is not None else np.nan),
            )

            # advance
            s = flatten_obs(next_obs)
            prev_loads = np.asarray(getattr(next_obs, "load_p")).copy()
            if done:
                break

        # ---- Episode-level aggregates
        def _pct(a, q):
            return float(np.percentile(a, q)) if len(a) else float("nan")

        ep_meta = {
            "system": pathlib.Path(grid_xml).stem,
            "agent": "DQN-deployed",
            "episode_seconds": float(getattr(next_obs, "current_step", 0.0)),

            # Forward-pass latency ℓ_act (ms)
            "act_ms_mean": float(np.mean(infer_ms_acc)) if infer_ms_acc else float("nan"),
            "act_ms_p50": _pct(infer_ms_acc, 50),
            "act_ms_p95": _pct(infer_ms_acc, 95),

            # Closed-loop latency ℓ_loop (ms)
            "loop_ms_mean": float(np.mean(loop_ms_acc)) if loop_ms_acc else float("nan"),
            "loop_ms_p50": _pct(loop_ms_acc, 50),
            "loop_ms_p95": _pct(loop_ms_acc, 95),

            # Env step (ms) for reference
            "env_ms_mean": float(np.mean(env_ms_acc)) if env_ms_acc else float("nan"),

            # Real-time factor ρ
            "rho_mean": float(np.mean(rho_vals)) if rho_vals else float("nan"),
            "rho_p50": _pct(rho_vals, 50),
            "rho_p05": _pct(rho_vals, 5),

            # Memory footprint
            "policy_size_mib": float(model_mib),
            "proc_rss_peak_mib": float(max(rss_mib_series)) if rss_mib_series else float("nan"),
        }

        if device.type == "cuda":
            ep_meta["cuda_peak_alloc_mib"] = float(torch.cuda.max_memory_allocated() / (1024.0 * 1024.0))
            ep_meta["cuda_peak_reserved_mib"] = float(torch.cuda.max_memory_reserved() / (1024.0 * 1024.0))
            torch.cuda.reset_peak_memory_stats()

        timelog.end_episode(meta=ep_meta)

    # -------- Teardown --------
    timelog.start("logging_s")
    saver.close()
    env.close()
    timelog.stop("logging_s")

    # Episode summaries + section totals
    timelog.dump_summaries()
    print(f"[deploy] done. Logs in {pathlib.Path(outdir) / 'deploy'}")

def main():
    parser = argparse.ArgumentParser(description="GridPACK-Grid2Op RL pipeline")
    parser.add_argument("--grid_xml", choices=["input_9bus.xml","input_9bus_fault.xml","input_39bus_IBR.xml","input_39bus_IBR_fault.xml","input_145.xml","input_145_fault.xml"], default="input_9bus.xml",
                        help="input_9bus.xml,input_9bus_fault.xml,input_39bus_IBR.xml,input_39bus_IBR_fault.xml,input_145.xml,input_145_fault.xml")
    parser.add_argument("--mode", choices=["simulate","train","deploy"], default="simulate",
                        help="simulate: do-nothing; train: learn DQN; deploy: run trained DQN")
    parser.add_argument("--outdir", type=str, default="runs/dqn_grid2op", help="Output directory")
    parser.add_argument("--episodes", type=int, default=1, help="Episodes for simulate/deploy")
    parser.add_argument("--max-steps", type=int, default=5, help="Max steps per episode for simulate/deploy")
    parser.add_argument("--total-steps", type=int, default=TOTAL_STEPS, help="Total env steps for training")
    parser.add_argument("--model-path", type=str, default="./models/dqn_grid2op.pt", help="Path to trained model for deploy")
    parser.add_argument("--log-sim-during-train", action="store_true", help="Also log per-step arrays during training")
    parser.add_argument("--seed", type=int, default=123, help="Random seed")
    args = parser.parse_args()

    grid_xml = args.grid_xml
    config_path = f"{CONFIG_DIR}/grid2op_{grid_xml.split('.')[0]}"

    # Repro
    random.seed(args.seed); np.random.seed(args.seed)
    torch.manual_seed(args.seed)

    print(f"Grid2Op {grid2op.__version__} | Torch {torch.__version__} | Mode={args.mode}")

    outdir = os.path.join(args.outdir, os.path.basename(grid_xml).split(".")[0])
    os.makedirs(outdir, exist_ok=True)

    if args.mode == "simulate":
        run_simulation(grid_xml=grid_xml, config_path=config_path, outdir=outdir, episodes=args.episodes, max_steps=args.max_steps, log_each_step=True)
    elif args.mode == "train":
        train(grid_xml=grid_xml, config_path=config_path, outdir=outdir, max_steps=args.max_steps, total_steps=args.total_steps, log_each_step=args.log_sim_during_train)
    elif args.mode == "deploy":
        deploy(grid_xml=grid_xml, config_path=config_path, model_path=args.model_path, outdir=outdir, episodes=args.episodes, max_steps=args.max_steps)

if __name__ == "__main__":
    main()