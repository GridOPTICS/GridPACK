import os, sys, time, csv, json, pathlib
from typing import Optional, Dict, List, Union, Any
from collections import deque
import psutil
import statistics

import numpy as np
import pandas as pd
from config import T_FINAL_CHECK, T_FAULT_CLEAR

def _proc():
    return psutil.Process(os.getpid())

def _rss_mib() -> float:
    return _proc().memory_info().rss / (1024 * 1024)

class TimingLogger:
    """
    Lightweight timing + memory logger for:
      - section totals (setup/training/inference/logging/simulate)
      - per-step env latency + RSS
      - episode-level summaries (avg, p95 latency; peak & delta RSS)
    """
    def __init__(self, outdir: str, run_name: str, enable_traces: bool = True):
        self.outdir = pathlib.Path(outdir) / run_name
        self.outdir.mkdir(parents=True, exist_ok=True)
        self.enable_traces = enable_traces

        # sections
        self.sections: Dict[str, float] = {
            "setup_s": 0.0,
            "training_s": 0.0,
            "inference_s": 0.0,
            "logging_s": 0.0,
            "simulate_s": 0.0,
        }
        self._t0: Optional[float] = None
        self._active: Optional[str] = None

        # per-episode / per-step
        self.current_episode: Optional[int] = None
        self.step_lat_ms: List[float] = []
        self.step_rss_mib: List[float] = []
        self._step_rows: List[dict] = []

        # summaries across episodes
        self.episodes: List[Dict[str, Any]] = []

    # ---- sections
    def start(self, section: str):
        if self._active is not None:
            self.stop(self._active)
        self._active = section
        self._t0 = time.perf_counter()

    def stop(self, section: str):
        if self._active != section or self._t0 is None:
            return
        dt = time.perf_counter() - self._t0
        self.sections[section] = self.sections.get(section, 0.0) + dt
        self._active, self._t0 = None, None

    # ---- episodes
    def begin_episode(self, ep_idx: int):
        self.current_episode = ep_idx
        self.step_lat_ms.clear()
        self.step_rss_mib.clear()
        self._step_rows.clear()
        # capture baseline memory at episode start
        self._rss0 = _rss_mib()

    def add_step(self, latency_ms: float, rss_mib: float, **kwargs):
        self.step_lat_ms.append(latency_ms)
        self.step_rss_mib.append(rss_mib)

        row = {
            "step": len(self._step_rows),
            # keep legacy names so old readers still work
            "latency_ms": float(latency_ms),
            "rss_mib": float(rss_mib),
        }
        # merge new fields
        for k, v in kwargs.items():
            # make floats where sensible
            row[k] = float(v) if isinstance(v, (int, float, np.floating)) and not isinstance(v, bool) else v
        self._step_rows.append(row)

    def end_episode(self, meta: Optional[Dict[str, Any]] = None):
        if not self.step_lat_ms:
            return
        lat_mean = float(sum(self.step_lat_ms) / len(self.step_lat_ms))
        lat_p95 = float(statistics.quantiles(self.step_lat_ms, n=100)[94]) if len(self.step_lat_ms) >= 20 else max(self.step_lat_ms)
        peak_rss = float(max(self.step_rss_mib)) if self.step_rss_mib else _rss_mib()
        delta_rss = float(peak_rss - getattr(self, "_rss0", peak_rss))

        summary = {
            "episode": self.current_episode,
            "n_steps": len(self.step_lat_ms),
            "latency_ms_avg": lat_mean,
            "latency_ms_p95": lat_p95,
            "peak_rss_mib": peak_rss,
            "delta_rss_mib": delta_rss,
        }
        if meta:
            summary.update(meta)
        self.episodes.append(summary)

        if self.enable_traces:
            tracef = self.outdir / f"steps_ep{self.current_episode:05d}.csv"
            all_keys = []
            seen = set()
            for r in self._step_rows:
                for k in r.keys():
                    if k not in seen:
                        seen.add(k); all_keys.append(k)
            with tracef.open("w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=all_keys)
                w.writeheader()
                for r in self._step_rows:
                    w.writerow(r)

        # reset
        self.current_episode = None
        self.step_lat_ms.clear()
        self.step_rss_mib.clear()
        self._step_rows.clear()

    # ---- dump summaries
    def dump_summaries(self):
        # per-episode summaries
        epf = self.outdir / "system_efficiency_episodes.csv"
        with epf.open("w", newline="") as f:
            if self.episodes:
                keys = list(self.episodes[0].keys())
            else:
                keys = ["episode","n_steps","latency_ms_avg","latency_ms_p95","peak_rss_mib","delta_rss_mib"]
            w = csv.DictWriter(f, fieldnames=keys)
            w.writeheader()
            for row in self.episodes:
                w.writerow(row)

        # section totals
        totf = self.outdir / "section_totals.json"
        with totf.open("w") as f:
            json.dump(self.sections, f, indent=2)

              
# ----------------- Logging -----------------
class TrainLogger:
    def __init__(self, outdir="runs/dqn"):
        self.outdir = pathlib.Path(outdir)
        (self.outdir / "ckpts").mkdir(parents=True, exist_ok=True)
        # rolling
        self.losses_recent = deque(maxlen=100)
        # files
        self.step_csv = open(self.outdir / "steps.csv", "w", newline="")
        self.ep_csv = open(self.outdir / "episodes.csv", "w", newline="")
        self.step_w = csv.DictWriter(self.step_csv,
            fieldnames=["global_step","episode","t","reward","action","epsilon","loss"])
        self.ep_w = csv.DictWriter(self.ep_csv,
            fieldnames=["episode","steps_total","return","len","epsilon_end","avg_loss100","wall_sec"])
        self.step_w.writeheader()
        self.ep_w.writeheader()
        self.start_time = time.time()

    def log_step(self, global_step, episode, t, reward, action, epsilon, loss):
        if loss is not None:
            self.losses_recent.append(float(loss))
        self.step_w.writerow({
            "global_step": global_step,
            "episode": episode,
            "t": t,
            "reward": float(reward),
            "action": int(action),
            "epsilon": float(epsilon),
            "loss": (None if loss is None else float(loss)),
        })

    def log_episode(self, episode, steps_total, ep_return, ep_len, epsilon_end):
        wall = time.time() - self.start_time
        self.ep_w.writerow({
            "episode": episode,
            "steps_total": steps_total,
            "return": float(ep_return),
            "len": ep_len,
            "epsilon_end": float(epsilon_end),
            "avg_loss100": (None if not self.losses_recent else float(np.mean(self.losses_recent))),
            "wall_sec": float(wall),
        })

    def flush(self):
        self.step_csv.flush()
        self.ep_csv.flush()

    def close(self):
        self.step_csv.close()
        self.ep_csv.close()

def _flat_col(prefix: str, id_levels: Union[str, int, tuple, list], keep_col: str) -> str:
    """
    Build "<prefix><id_levels>_<keep_col>" with whitespace stripped.
    - id_levels can be a scalar (e.g., '5') or a tuple/list (e.g., (5, 'G1'))
    """
    if isinstance(id_levels, (list, tuple)):
        id_part = "_".join(str(x).strip() for x in id_levels)
    else:
        id_part = str(id_levels).strip()
    return f"{prefix}{id_part}_{str(keep_col).strip()}"


class SimDataSaver:
    """
    Real-time single-file logger.
    - On __init__, freezes IDs and id_col for gen/load/bus.
    - On each log_step(), pivots res_* to wide form, merges, adds scalars/arrays, writes.
    - Columns: scalars + arrays + wide pivots (fixed).
    """

    def __init__(self, g, outdir: str, run_name: str,
                 n_bus: Optional[int] = None, n_load: Optional[int] = None,
                 attrs_present: Optional[Dict[str, bool]] = None):

        base = pathlib.Path(outdir) / run_name
        base.mkdir(parents=True, exist_ok=True)
        self.base = base
        self.attrs_present = attrs_present or {}

        # ---- Freeze ID_COLS ----
        self.gen_id_col = ["bus", "name"]
        self.load_id_col = ["bus", "name"]
        self.bus_id_col = "bus_id"

        # breakpoint()
        # ---- Freeze IDs (unique values) ----
        self.gen_id_pairs = []
        if hasattr(g, "gen_logger") and g.gen_logger:
            gen_frame = pd.concat(g.gen_logger, ignore_index=True)
            if set(self.gen_id_col) <= set(gen_frame.columns):
                self.gen_id_pairs = (
                    gen_frame[self.gen_id_col]
                    .drop_duplicates()
                    .apply(lambda r: (getattr(g, "BUS_MAPPING_LOCAL2REAL", {}).get(r["bus"], r["bus"]),
                                      r["name"]), axis=1).tolist()
                )
        self.gen_params = ["name", "bus", "p_mw", "q_mvar", "vn_kv"]

        self.load_id_pairs = []
        if hasattr(g, "load_logger") and g.load_logger:
            load_frame = pd.concat(g.load_logger, ignore_index=True)
            if set(self.load_id_col) <= set(load_frame.columns):
                self.load_id_pairs = (
                    load_frame[self.load_id_col]
                    .drop_duplicates()
                    .apply(lambda r: (getattr(g, "BUS_MAPPING_LOCAL2REAL", {}).get(r["bus"], r["bus"]),
                                      r["name"]), axis=1).tolist()
                )
        # breakpoint()
        self.load_params = ["name", "bus", "p_mw", "q_mvar"]

        self.bus_ids = []
        if hasattr(g, "bus_logger") and g.bus_logger:
            bus_frame = pd.concat(g.bus_logger, ignore_index=True)
            if "id" in bus_frame.columns:
                self.bus_ids = list(pd.unique(bus_frame["id"]))
        self.bus_params = ["base_kv", "case_sbase", "vn_pu", "vn_kv", "frequency"]

        # ---- Scalar columns per step ----
        self._scalar_fields = [
            "episode", "t", "tick", "step", "action_idx", "reward",
            "v_threshold", "v_pen", "shed_pen", "invalid_pen",
            "t_since_fault", "fault_cleared"
        ]

        # ---- Array columns (obs vectors) ----
        self._sizes: Dict[str, int] = {}
        if n_bus is not None and self.attrs_present.get("bus_v", False):
            self._sizes["bus_v"] = int(n_bus)
        if n_bus is not None and self.attrs_present.get("bus_freq", False):
            self._sizes["bus_freq"] = int(n_bus)
        if n_load is not None and self.attrs_present.get("load_v", False):
            self._sizes["load_v"] = int(n_load)
        if n_load is not None and self.attrs_present.get("load_p", False):
            self._sizes["load_p"] = int(n_load)
        if n_load is not None and self.attrs_present.get("load_q", False):
            self._sizes["load_q"] = int(n_load)

        array_cols = [f"{name}_{i}" for name, sz in self._sizes.items() for i in range(sz)]

        # ---- Wide columns (fixed once, consistent with _to_map) ----
        self._wide_cols = []
        for pair in self.gen_id_pairs:
            for p in self.gen_params:
                self._wide_cols.append(_flat_col("gen_", pair, p))
        for pair in self.load_id_pairs:
            for p in self.load_params:
                self._wide_cols.append(_flat_col("load_", pair, p))
        for bid in self.bus_ids:
            for p in self.bus_params:
                self._wide_cols.append(_flat_col("bus_", bid, p))

        # ---- Final header ----
        self._header = self._scalar_fields + array_cols + self._wide_cols

        # ---- Open CSV ----
        self._path = self.base / "sim_all.csv"
        self._file = open(self._path, "w", newline="")
        self._writer = csv.DictWriter(self._file, fieldnames=self._header)
        self._writer.writeheader()

        # Buffer latest arrays until log_step()
        self._last_arrays: Dict[str, Optional[np.ndarray]] = {k: None for k in self._sizes.keys()}

    # ---- Helpers ----
    def _arrsafe(self, obs, name: str) -> Optional[np.ndarray]:
        if not self.attrs_present.get(name, False): return None
        if not hasattr(obs, name): return None
        arr = getattr(obs, name)
        if arr is None: return None
        a = np.asarray(arr).ravel().astype(float, copy=False)
        sz = self._sizes.get(name, None)
        if sz is not None:
            if a.size > sz: a = a[:sz]
            elif a.size < sz:
                pad = np.full(sz - a.size, np.nan, dtype=float)
                a = np.concatenate([a, pad])
        return a

    def _to_map(self, df: pd.DataFrame, id_col, keep_cols: List[str], prefix: str,
                tick_col: str = "tick", bus_map=None) -> pd.DataFrame:
        """
        Return a WIDE DF: index=tick, cols="<prefix><id_levels>_<keep_col>"
        """
        if df is None or len(df) == 0:
            return pd.DataFrame()

        d = df.copy()
        if bus_map is not None and isinstance(id_col, (list, tuple)) and "bus" in d.columns:
            d["bus"] = d["bus"].map(lambda x: bus_map.get(x, x))

        d[tick_col] = pd.to_numeric(d[tick_col], errors="coerce")
        d = d.dropna(subset=[tick_col])

        wide = d.pivot(index=tick_col, columns=id_col, values=keep_cols).sort_index()

        # flatten
        new_cols = []
        for tup in wide.columns.to_flat_index():
            keep = str(tup[0]).strip()
            ids = tup[1:]
            new_cols.append(_flat_col(prefix, ids if len(ids) != 1 else ids[0], keep))
        wide.columns = new_cols
        wide.index.name = "tick"
        return wide

    # ---- Step I/O ----
    def log_arrays(self, obs):
        for name in self._sizes.keys():
            self._last_arrays[name] = self._arrsafe(obs, name)

    def log_step(self, g, episode: int, t: float, step: int, a_idx: int, reward: float, rc: Dict[str, float]):
        gen_wide = load_wide = bus_wide = pd.DataFrame()

        if getattr(g, "gen_logger", None):
            gen_frame = pd.concat(g.gen_logger, ignore_index=True)
            gen_wide = self._to_map(gen_frame, id_col=self.gen_id_col, keep_cols=self.gen_params, prefix="gen_", bus_map=getattr(g, "BUS_MAPPING_LOCAL2REAL", None))
        if getattr(g, "load_logger", None):
            load_frame = pd.concat(g.load_logger, ignore_index=True)
            load_wide = self._to_map(load_frame, id_col=self.load_id_col, keep_cols=self.load_params, prefix="load_", bus_map=getattr(g, "BUS_MAPPING_LOCAL2REAL", None))
        if getattr(g, "bus_logger", None):
            bus_frame = pd.concat(g.bus_logger, ignore_index=True).rename(columns={"id": "bus_id"})
            bus_wide = self._to_map(bus_frame, id_col=self.bus_id_col, keep_cols=self.bus_params, prefix="bus_")

        # Merge on tick
        parts = [df for df in (gen_wide, load_wide, bus_wide) if not df.empty]
        if not parts:
            return
        merged = parts[0].join(parts[1:], how="outer") if len(parts) > 1 else parts[0]
        merged = merged.sort_index().reset_index()

        # Scalars
        scalars = {
            "episode": int(episode),
            "t": float(t),
            "step": int(step),
            "action_idx": (None if a_idx is None else int(a_idx)),
            "reward": float(reward),
            "v_threshold": float(rc.get("threshold", np.nan)),
            "v_pen": float(rc.get("v_pen", np.nan)),
            "shed_pen": float(rc.get("shed_pen", np.nan)),
            "invalid_pen": float(rc.get("invalid_pen", np.nan)),
            "t_since_fault": float(rc.get("t_since_fault", np.nan)),
            "fault_cleared": int(1 if rc.get("fault_cleared", False) else 0),
        }
        merged = merged.assign(**scalars)

        # Arrays
        array_cols = [f"{name}_{i}" for name, sz in self._sizes.items() for i in range(sz)]
    
        if array_cols:
            arr_block = pd.DataFrame(np.nan, index=merged.index, columns=array_cols)

            # Fill only the row matching current tick (vectorized row assignment)
            curr_tick = float(rc.get("tick", np.nan))
            # ensure numeric comparison; if tick column might be str, coerce:
            tick_vals = pd.to_numeric(merged["tick"], errors="coerce")
            sel = np.isclose(tick_vals.values.astype(float), curr_tick, equal_nan=False)

            if sel.any():
                row_idx = np.flatnonzero(sel)[0]
                # Assign each array slice to its columns in that single row
                for name, sz in self._sizes.items():
                    arr = self._last_arrays.get(name)
                    if arr is not None:
                        cols = [f"{name}_{i}" for i in range(sz)]
                        # assign all sz values at once
                        arr_block.loc[row_idx, cols] = np.asarray(arr, dtype=float)[:sz]

            # One concat -> no fragmentation
            merged = pd.concat([merged, arr_block], axis=1)

        # Reorder & write
        merged = merged.reindex(columns=self._header, fill_value=np.nan)
        self._writer.writerows(merged.to_dict(orient="records"))

    def flush(self):
        try:
            self._file.flush()
        except Exception:
            pass

    def close(self):
        try:
            self._file.flush()
            self._file.close()
        except Exception:
            pass

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