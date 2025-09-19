import os, sys, math, random
sys.path += [
    "/qfs/projects/gridpack_wind/grid2op_interface/grid2op_local/",
    "/qfs/projects/gridpack_wind/grid2op_interface/GridPACK/python/src/",
]

import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from datetime import timedelta

import pathlib
import argparse

from config import CONFIG_PATH, GRID_XML, TIME_STEP_S, LR, GAMMA, TARGET_SYNC, BATCH_SIZE, EPS_START, EPS_END, EPS_DECAY, CONTROLLED_LOADS, SHED_RATIO, TOTAL_STEPS

from model import MLP, Replay, build_action_list

from data_loggers import TrainLogger, SimDataSaver
from reward import FIDVRReward, reward_components
from utils import read_gridpack_stepsize, flatten_obs, get_sizes, detect_attrs



import grid2op
from grid2op import make
from grid2op.Action import BaseAction
from grid2op.Reward import BaseReward
from grid2op.Runner import Runner

from grid2op_backend import GridPACKBackend

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

def run_simulation(outdir: str, episodes: int, max_steps: int, log_each_step: bool=True):
    """Do-nothing agent (no shedding), just run and log."""
    env = make_env()
    # Prepare logging based on first obs
    obs0 = env.reset()
    attrs_present = detect_attrs(obs0)
    n_bus, n_load = get_sizes(obs0)
    saver = SimDataSaver(
        env.backend, 
        outdir, 
        run_name="simulate", 
        n_bus=n_bus, 
        n_load=n_load, 
        attrs_present=attrs_present
    )

    # Build an action list but use index 0 as "do nothing"
    actions = build_action_list(env, CONTROLLED_LOADS[GRID_XML], SHED_RATIO)
    do_nothing_idx = 0

    for ep in range(1, episodes+1):
        obs = env.reset()
        # init for reward components
        init_loads = np.asarray(getattr(obs, "load_p")).copy()
        prev_loads = init_loads.copy()

        for t in range(max_steps):
            a_idx = do_nothing_idx
            next_obs, r, done, info = env.step(actions[a_idx])

            # log arrays
            if log_each_step:
                rc = reward_components(next_obs, prev_loads, init_loads)
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

            prev_loads = np.asarray(getattr(next_obs, "load_p")).copy()
            if done: break

    saver.close()
    env.close()
    print(f"[simulate] done. Logs in {pathlib.Path(outdir) / 'simulate'}")

def train(outdir: str, total_steps: int, log_each_step: bool=False):
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    env = make_env()
    action_list = build_action_list(env, CONTROLLED_LOADS[GRID_XML], SHED_RATIO)

    o0 = env.reset()
    s0 = flatten_obs(o0)
    in_dim, out_dim = s0.size, len(action_list)

    q = MLP(in_dim, out_dim).to(device)
    qt = MLP(in_dim, out_dim).to(device)
    qt.load_state_dict(q.state_dict())
    opt = optim.Adam(q.parameters(), lr=LR)
    replay = Replay(100000)

    logger = TrainLogger(outdir=str(pathlib.Path(outdir) / "train"))
    SAVE_EVERY_STEPS = 1000

    # Optional sim saver during training
    saver = None
    if log_each_step:
        attrs_present = detect_attrs(o0)
        n_bus, n_load = get_sizes(o0)
        saver = SimDataSaver(outdir, run_name="train_sim", n_bus=n_bus, n_load=n_load, attrs_present=attrs_present)

    global_step = 0
    ep_idx = 0
    best_return = -1e9

    def epsilon(step):
        # linear decay
        if step >= EPS_DECAY: return EPS_END
        return EPS_START + (EPS_END - EPS_START) * (step / EPS_DECAY)

    while global_step < total_steps:
        ep_idx += 1
        obs = env.reset()
        s = flatten_obs(obs)
        ep_ret = 0.0
        done = False
        max_steps = int(math.ceil(300.0 / TIME_STEP_S))
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

            next_obs, r, d, info = env.step(action_list[a_idx])
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

            logger.log_step(global_step, ep_idx, t, r, a_idx, eps, (loss_value if loss_value is not None else None))

            if saver is not None:
                rc = reward_components(next_obs, prev_loads, init_loads)
                saver.log_arrays(next_obs)
                saver.log_step(ep_idx, next_obs.current_time, t, a_idx, float(r), rc)

            # target sync
            if global_step % TARGET_SYNC == 0:
                qt.load_state_dict(q.state_dict())

            # periodic checkpoint + flush
            if global_step % SAVE_EVERY_STEPS == 0:
                torch.save(q.state_dict(), pathlib.Path(outdir) / "train" / f"ckpts/step_{global_step}.pt")
                logger.flush()

            prev_loads = np.asarray(getattr(next_obs, "load_p")).copy()

            if d or global_step >= total_steps:
                print(f"ep {ep_idx:04d} | steps {global_step} | return {ep_ret:.3f} | eps {eps:.3f}")
                if ep_ret > best_return:
                    best_return = ep_ret
                    os.makedirs("./models", exist_ok=True)
                    torch.save(q.state_dict(), "./models/dqn_grid2op.pt")
                logger.log_episode(ep_idx, global_step, ep_ret, ep_len, last_eps)
                logger.flush()
                break

    env.close()
    if saver is not None: saver.close()
    logger.close()
    print("[train] done.")

def deploy(model_path: str, outdir: str, episodes: int, max_steps: int):
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    env = make_env()
    actions = build_action_list(env, CONTROLLED_LOADS[GRID_XML], SHED_RATIO)

    # Model dims from env
    o0 = env.reset()
    in_dim, out_dim = flatten_obs(o0).size, len(actions)
    q = MLP(in_dim, out_dim).to(device)
    q.load_state_dict(torch.load(model_path, map_location=device))
    q.eval()

    attrs_present = detect_attrs(o0)
    n_bus, n_load = get_sizes(o0)
    saver = SimDataSaver(outdir, run_name="deploy", n_bus=n_bus, n_load=n_load, attrs_present=attrs_present)

    for ep in range(1, episodes+1):
        obs = env.reset()
        s = flatten_obs(obs)
        init_loads = np.asarray(getattr(obs, "load_p")).copy()
        prev_loads = init_loads.copy()

        for t in range(max_steps):
            with torch.no_grad():
                qs = q(torch.from_numpy(s).unsqueeze(0).to(device))
                a_idx = int(torch.argmax(qs, dim=1).item())

            next_obs, r, done, info = env.step(actions[a_idx])
            s = flatten_obs(next_obs)

            rc = reward_components(next_obs, prev_loads, init_loads)
            saver.log_arrays(next_obs)
            saver.log_step(ep, next_obs.current_time, t, a_idx, float(r), rc)

            prev_loads = np.asarray(getattr(next_obs, "load_p")).copy()
            if done: break

    saver.close()
    env.close()
    print(f"[deploy] done. Logs in {pathlib.Path(outdir) / 'deploy'}")

def main():
    parser = argparse.ArgumentParser(description="GridPACK-Grid2Op RL pipeline")
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

    # Repro
    random.seed(args.seed); np.random.seed(args.seed)
    torch.manual_seed(args.seed)

    print(f"Grid2Op {grid2op.__version__} | Torch {torch.__version__} | Mode={args.mode}")

    outdir = os.path.join(args.outdir, os.path.basename(GRID_XML).split(".")[0])
    os.makedirs(outdir, exist_ok=True)

    if args.mode == "simulate":
        run_simulation(outdir=outdir, episodes=args.episodes, max_steps=args.max_steps, log_each_step=True)
    elif args.mode == "train":
        train(outdir=outdir, total_steps=args.total_steps, log_each_step=args.log_sim_during_train)
    elif args.mode == "deploy":
        deploy(model_path=args.model_path, outdir=outdir, episodes=args.episodes, max_steps=args.max_steps)

if __name__ == "__main__":
    main()