#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os,re, argparse, json
from pathlib import Path
from typing import Optional, List, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ---------- Global plot style ----------
plt.rcParams.update({
    "font.size": 14,
    "axes.titlesize": 16,
    "axes.labelsize": 14,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
})

# ---------- IO helpers ----------
def read_csv(p: Path) -> Optional[pd.DataFrame]:
    try:
        if p.exists():
            return pd.read_csv(p)
    except Exception as e:
        print(f"[warn] Could not read {p}: {e}")
    return None

def read_json(p: Path) -> Optional[dict]:
    try:
        if p.exists():
            return json.loads(p.read_text())
    except Exception as e:
        print(f"[warn] Could not read {p}: {e}")
    return None

# --------- helpers to autodetect columns ---------
def _cols(df, suffix):
    return [c for c in df.columns if c.endswith(suffix)]

def bus_voltage_cols(df):     return _cols(df, "_vn_pu")
def bus_frequency_cols(df):   return _cols(df, "_frequency")
def load_p_cols(df):          return [c for c in df.columns if c.startswith("load_") and c.endswith("_p_mw")]

# ---------- Inference (deploy) ----------
def summarize_deploy_inference(df: pd.DataFrame, label: str) -> dict:
    """
    Expect columns created in deploy ep_meta:
      act_ms_mean, act_ms_p50, act_ms_p95,
      loop_ms_mean, loop_ms_p50, loop_ms_p95,
      env_ms_mean,
      rho_mean, rho_p50, rho_p05,
      policy_size_mib, proc_rss_peak_mib,
      (optional) cuda_peak_alloc_mib, cuda_peak_reserved_mib
    """
    def _mstd(col):
        if col not in df.columns:
            return float("nan"), 0.0
        s = df[col].astype(float)
        m = float(s.mean()) if len(s) else float("nan")
        v = float(s.std(ddof=1)) if len(s) > 1 else 0.0
        return m, v

    out = {"label": label}
    keys = [
        "act_ms_mean", "act_ms_p50", "act_ms_p95",
        "loop_ms_mean", "loop_ms_p50", "loop_ms_p95",
        "env_ms_mean",
        "rho_mean", "rho_p50", "rho_p05",
        "policy_size_mib", "proc_rss_peak_mib",
        "cuda_peak_alloc_mib", "cuda_peak_reserved_mib",
    ]
    for k in keys:
        out[f"{k}_mean"], out[f"{k}_std"] = _mstd(k)
    # A convenience: episode wall vs sim time ratio (ρ_ep) if available
    if {"episode_seconds"}.issubset(df.columns) and "loop_ms_mean" in df.columns:
        # Approx wall time per episode ≈ n_steps * loop_ms_mean; we don’t have n_steps per ep here,
        # so just report mean(loop_ms_mean) alongside rho_mean already computed per-step in deploy.
        pass
    return out

def latex_table_deploy_inference(rows: List[dict]) -> str:
    """
    Build a LaTeX table for deploy (inference) metrics.

    - Matches the example layout when there's exactly one system:
        * specific row ordering
        * no CUDA rows
        * tabular {lc}
    - For multiple systems, expands to extra columns automatically.
    """
    if not rows:
        return ""

    systems = [r.get("label", f"System {i+1}") for i, r in enumerate(rows)]

    def fmt(x):
        try:
            if x is None or (isinstance(x, float) and np.isnan(x)):
                return r"--"
            return f"{float(x):.2f}"
        except Exception:
            return r"--"

    # Row spec: (title, mean_key_in_rows, std_key_in_rows)
    # Order matches your example.
    metrics = [
        (r"Forward-pass latency $\ell_{\text{act}}$ (ms)", "act_ms_mean_mean", "act_ms_mean_std"),
        (r"$\ell_{\text{act}}$ p95 (ms)",                  "act_ms_p95_mean",  "act_ms_p95_std"),
        (r"Closed-loop latency $\ell_{\text{loop}}$ (ms)", "loop_ms_mean_mean","loop_ms_mean_std"),
        (r"$\ell_{\text{loop}}$ p95 (ms)",                 "loop_ms_p95_mean", "loop_ms_p95_std"),
        (r"Real-time factor $\rho$ (—)",                   "rho_mean_mean",    "rho_mean_std"),
        (r"$\rho$ p50 (—)",                                "rho_p50_mean",     "rho_p50_std"),
        (r"$\rho$ p05 (—)",                                "rho_p05_mean",     "rho_p05_std"),
        (r"Env step latency (ms)",                         "env_ms_mean_mean", "env_ms_mean_std"),
        (r"Policy size (MiB)",                             "policy_size_mib_mean","policy_size_mib_std"),
        (r"Peak RSS (MiB)",                                "proc_rss_peak_mib_mean","proc_rss_peak_mib_std"),
    ]

    # Header
    header = "\\begin{table}[t]\n\\renewcommand{\\arraystretch}{1.3}\n\\centering\n"
    header += "\\caption{Inference performance (deploy): latency, real-time factor, and memory (mean $\\pm$ std across episodes).}\n"
    header += "\\label{tab:deploy_inference}\n"

    # Column spec: {lc} for single system, otherwise {lccc...}
    colspec = "l" + "c" * len(systems)
    if len(systems) == 1:
        colspec = "lc"
    header += f"\\begin{{tabular}}{{{colspec}}}\n\\hline\n"
    header += "\\textbf{Metric} & " + " & ".join(systems) + " \\\\\n\\hline\n"

    # Body
    body_lines = []
    for title, mean_key, std_key in metrics:
        cells = []
        for r in rows:
            m = fmt(r.get(mean_key))
            s = fmt(r.get(std_key))
            cells.append(f"{m} $\\pm$ {s}")
        body_lines.append(f"{title} & " + " & ".join(cells) + r" \\")
    body = "\n".join(body_lines) + "\n"

    # Footer
    footer = "\\hline\n\\end{tabular}\n\\end{table}\n"
    return header + body + footer

def plot_deploy_lines(ep_df: pd.DataFrame, outdir: Path, suffix: str):
    """
    Simple per-episode line plots: act_ms_mean, loop_ms_mean, rho_mean, and memory peaks.
    """
    outdir.mkdir(parents=True, exist_ok=True)
    if "episode" not in ep_df.columns:
        # try to synthesize episode index
        ep_df = ep_df.reset_index().rename(columns={"index":"episode"})
    df = ep_df.sort_values("episode")

    def line(col, ylabel, fname):
        if col not in df.columns: 
            print(f"[info] deploy: missing column '{col}', skipping {fname}")
            return
        plt.figure()
        plt.plot(df["episode"], df[col], marker="o", linewidth=1)
        plt.xlabel("Episode")
        plt.ylabel(ylabel)
        plt.tight_layout()
        plt.savefig(outdir / f"{fname}_{suffix}.png", dpi=200)
        plt.close()

    line("act_ms_mean", "Forward-pass latency (ms)", "deploy_act_ms_mean")
    line("loop_ms_mean", "Closed-loop latency (ms)", "deploy_loop_ms_mean")
    line("rho_mean", "Real-time factor (—)", "deploy_rho_mean")

    # Memory (bars)
    for mem_col, ylabel, fname in [
        ("proc_rss_peak_mib", "Peak RSS (MiB)", "deploy_peak_rss"),
        ("cuda_peak_alloc_mib", "CUDA peak alloc (MiB)", "deploy_cuda_peak_alloc"),
        ("cuda_peak_reserved_mib", "CUDA peak reserved (MiB)", "deploy_cuda_peak_reserved"),
    ]:
        if mem_col not in df.columns: 
            continue
        plt.figure()
        plt.bar(df["episode"], df[mem_col])
        plt.xlabel("Episode")
        plt.ylabel(ylabel)
        plt.tight_layout()
        plt.savefig(outdir / f"{fname}_{suffix}.png", dpi=200)
        plt.close()

def plot_deploy_step_hists_if_available(deploy_dir: Path, outdir: Path, suffix: str):
    """
    Optional: if step traces contain per-step 'infer_ms' or 'loop_ms', plot histograms.
    Falls back silently if not present.
    """
    files = sorted(deploy_dir.glob("steps_ep*.csv"))
    if not files:
        print(f"[info] deploy: no per-step traces found.")
        return

    def collect(col):
        arrs = []
        for f in files:
            df = read_csv(f)
            if df is None or df.empty or col not in df.columns:
                continue
            arrs.append(df[col].astype(float).values)
        return np.concatenate(arrs) if arrs else np.array([])

    for col, xlabel, fname in [
        ("infer_ms", "Forward-pass latency (ms)", "deploy_infer_ms_hist"),
        ("loop_ms",  "Closed-loop latency (ms)",  "deploy_loop_ms_hist"),
    ]:
        data = collect(col)
        if data.size == 0:
            print(f"[info] deploy: '{col}' not in step traces; skipping {fname}")
            continue
        plt.figure()
        plt.hist(data, bins=50)
        plt.xlabel(xlabel)
        plt.ylabel("Count")
        plt.tight_layout()
        plt.savefig(outdir / f"{fname}_{suffix}.png", dpi=200)
        plt.close()

# --------- core metrics ---------
def voltage_settling_time(df_ep, v_cols, thr_col="v_threshold", after_fault_only=True, consec=3):
    """
    First time (in seconds) after which ALL bus voltages stay >= v_threshold
    for 'consec' consecutive steps. Returns np.nan if never settles.
    Assumes 't' increases by 1s; otherwise dt is inferred from median diff.
    """
    if len(df_ep) == 0: return np.nan
    dt = float(np.median(np.diff(df_ep["t"].unique()))) if df_ep["t"].nunique() > 1 else 1.0
    # filter window: after clearing, if column available
    if after_fault_only and "fault_cleared" in df_ep.columns:
        # include rows at/after first clearing == 1
        idx_clear = df_ep.index[df_ep["fault_cleared"] == 1]
        if len(idx_clear):
            start_idx = idx_clear.min()
            df_ep = df_ep.loc[start_idx:]
    # threshold per-row (can vary); require all buses >= threshold each row
    ok = (df_ep[v_cols].min(axis=1) >= df_ep[thr_col])
    # find first index where we have 'consec' Trues in a row
    run = 0
    for i, val in enumerate(ok.values):
        run = (run + 1) if val else 0
        if run >= consec:
            # time at start of the satisfying run relative to episode start of this window
            t0 = df_ep["t"].iloc[i - consec + 1] - df_ep["t"].iloc[0]
            return float(t0)
    return np.nan

def load_shed_mwh(df_ep, load_cols):
    """
    Integrate (baseline MW - current MW)_+ over the episode -> MWh.
    Baseline = first row of the episode. dt inferred from 't' (s).
    """
    if len(df_ep) == 0: return np.nan
    dt_s = float(np.median(np.diff(df_ep["tick"].unique()))) if df_ep["tick"].nunique() > 1 else 1.0
    dt_h = dt_s / 3600.0
    base = df_ep[load_cols].iloc[0].astype(float).values
    cur  = df_ep[load_cols].astype(float).values
    shed_mw = np.maximum(0.0, (base - cur))  # shape: (T, nloads)
    # sum over loads then integrate over time
    total_mw = shed_mw.sum(axis=1)           # (T,)
    mwh = float((total_mw * dt_h).sum())
    return mwh

def evaluate_actions(df: pd.DataFrame, outdir: str | Path, no_op_idx: int = 0,
                     effect_window: int = 3, pre_window: int = 3):
    """
    Evaluate policy via chosen actions using sim_all.csv.
    Assumptions:
      - 'action_idx' == no_op_idx means 'do nothing'
      - 't' is step index in seconds (or constant step); 'fault_cleared' ∈ {0,1}
      - Bus voltages columns end with '_vn_pu'; bus frequencies end with '_frequency'
    """
    outdir = Path(outdir); outdir.mkdir(parents=True, exist_ok=True)
    # df = pd.read_csv(sim_all_csv)

    # ---- find voltage/frequency columns ----
    vcols = [c for c in df.columns if c.endswith("_vn_pu")]
    fcols = [c for c in df.columns if c.endswith("_frequency")]
    if not vcols or not fcols:
        raise ValueError("Expected *_vn_pu and *_frequency columns present.")

    # normalize sorting
    sort_cols = [c for c in ["episode", "t", "step"] if c in df.columns]
    df = df.sort_values(sort_cols)

    # per-episode metrics
    rows = []
    act_hist = {}
    for ep, dfe in df.groupby("episode", sort=True):
        dfe = dfe.copy()
        dfe["is_action"] = (dfe["action_idx"] != no_op_idx)

        # time to first action after fault clear
        t_clear_idx = dfe.index[dfe["fault_cleared"] == 1] if "fault_cleared" in dfe.columns else []
        t_clear = dfe.loc[t_clear_idx.min(), "t"] if len(t_clear_idx) else dfe["t"].min()
        post_clear = dfe[dfe["t"] >= t_clear]
        first_act_t = np.nan
        if not post_clear.empty:
            aidx = post_clear.index[post_clear["is_action"]]
            if len(aidx):
                first_act_t = float(post_clear.loc[aidx.min(), "t"] - t_clear)

        # actions per episode & inter-action gaps
        act_times = dfe.loc[dfe["is_action"], "t"].values
        n_actions = int(len(act_times))
        gaps = np.diff(act_times) if n_actions >= 2 else np.array([])

        # action distribution (episode-level)
        counts = dfe["action_idx"].value_counts()
        probs = (counts / counts.sum()).sort_index()
        entropy = float(-(probs * np.log(probs + 1e-12)).sum())

        # accumulate global histogram
        for k, v in counts.items():
            act_hist[k] = act_hist.get(k, 0) + int(v)

        # simple post-action effect estimate for last action in episode
        # (you could do this for every action and average)
        dv_mean = np.nan; dfreq_mean = np.nan
        if n_actions > 0:
            last_t = int(act_times[-1])
            pre = dfe[(dfe["t"] >= last_t - pre_window) & (dfe["t"] < last_t)]
            post = dfe[(dfe["t"] > last_t) & (dfe["t"] <= last_t + effect_window)]
            if not pre.empty and not post.empty:
                v_pre = pre[vcols].mean().mean()
                v_post = post[vcols].mean().mean()
                f_pre = pre[fcols].mean().mean()
                f_post = post[fcols].mean().mean()
                dv_mean = float(v_post - v_pre)
                dfreq_mean = float(f_post - f_pre)

        rows.append({
            "episode": ep,
            "first_action_delay_s": first_act_t,
            "n_actions": n_actions,
            "mean_gap_s": float(gaps.mean()) if gaps.size else np.nan,
            "action_entropy": entropy,
            "post_action_dV_mean": dv_mean,
            "post_action_dFreq_mean": dfreq_mean,
        })

    metrics = pd.DataFrame(rows).sort_values("episode")

    # ---- Save metrics CSV ----
    metrics_csv = outdir / "action_metrics_by_episode.csv"
    metrics.to_csv(metrics_csv, index=False)

    # ---- Global action histogram plot ----
    if act_hist:
        keys = sorted(act_hist.keys())
        vals = [act_hist[k] for k in keys]
        plt.figure()
        plt.bar(keys, vals)
        plt.xlabel("action_idx")
        plt.ylabel("count")
        plt.title("Action index distribution (all episodes)")
        plt.tight_layout()
        plt.savefig(outdir / "action_histogram.png", dpi=200)
        plt.close()

    # ---- Per-episode plots ----
    def lineplot(col, ylabel, fname):
        if col not in metrics.columns: return
        plt.figure()
        plt.plot(metrics["episode"].values, metrics[col].values, marker="o", linewidth=1)
        plt.xlabel("Episode")
        plt.ylabel(ylabel)
        plt.title(f"{ylabel} vs Episode")
        plt.tight_layout()
        plt.savefig(outdir / fname, dpi=200)
        plt.close()

    lineplot("first_action_delay_s", "First action delay (s)", "first_action_delay_vs_episode.png")
    lineplot("n_actions", "Actions per episode", "actions_per_episode.png")
    lineplot("mean_gap_s", "Mean inter-action gap (s)", "mean_interaction_gap_vs_episode.png")
    lineplot("action_entropy", "Action entropy", "action_entropy_vs_episode.png")
    lineplot("post_action_dV_mean", "Post-action ΔVoltage (p.u.)", "post_action_dV_vs_episode.png")
    lineplot("post_action_dFreq_mean", "Post-action ΔFrequency (Hz)", "post_action_dFreq_vs_episode.png")

    print(f"[done] Wrote {metrics_csv} and action plots to {outdir}")
    return metrics

def frequency_metrics(df_ep, f_cols, band=(59.5, 60.5)):
    """
    Return (nadir_hz, time_outside_band_s).
    Nadir = min across all buses & times in episode window.
    time_outside = seconds where any bus is outside [lo, hi].
    """
    if len(df_ep) == 0: return np.nan, np.nan
    dt_s = float(np.median(np.diff(df_ep["t"].unique()))) if df_ep["t"].nunique() > 1 else 1.0
    F = df_ep[f_cols].astype(float)
    nadir = float(F.min().min())
    lo, hi = band
    outside_mask = (F.lt(lo) | F.gt(hi)).any(axis=1)
    time_outside_s = float(outside_mask.sum() * dt_s)
    return nadir, time_outside_s

# --------- runner: per-episode + aggregate + LaTeX ---------
def summarize_policy_quality(df: pd.DataFrame):
    # df = pd.read_parquet(sim_all_pq)
    vcols = bus_voltage_cols(df)
    fcols = bus_frequency_cols(df)
    lcols = load_p_cols(df)

    if not len(vcols) or not len(fcols) or not len(lcols):
        raise ValueError("Missing expected columns: voltages(_vn_pu), frequency(_frequency), or loads(load_*_p_mw).")

    rows = []
    for ep, dfe in df.groupby("episode", sort=True):
        # breakpoint()
        dfe = dfe.sort_values("tick")
        v_settle_s = voltage_settling_time(dfe, vcols, thr_col="v_threshold", after_fault_only=True, consec=3)
        mwh = load_shed_mwh(dfe, lcols)
        nadir, t_out = frequency_metrics(dfe, fcols, band=(59.5, 60.5))
        rows.append({
            "episode": ep,
            "voltage_settle_s": v_settle_s,
            "load_shed_mwh": mwh,
            "freq_nadir_hz": nadir,
            "freq_time_outside_s": t_out,
        })
    res = pd.DataFrame(rows)
    # breakpoint()

    def m(x): return float(np.nanmean(x))
    def s(x): return float(np.nanstd(x, ddof=1)) if len(x.dropna())>1 else 0.0
    def p95(x): return float(np.nanpercentile(x, 95)) if len(x.dropna()) else np.nan

    agg = {
        "Voltage settling time to TVRC (s)": (m(res["voltage_settle_s"]), s(res["voltage_settle_s"]), p95(res["voltage_settle_s"])),
        "Load shed (MWh / episode)":         (m(res["load_shed_mwh"]), s(res["load_shed_mwh"]), p95(res["load_shed_mwh"])),
        "Frequency nadir (Hz)":              (m(res["freq_nadir_hz"]), s(res["freq_nadir_hz"]), p95(res["freq_nadir_hz"])),
        "Time outside [59.5,60.5] Hz (s)":   (m(res["freq_time_outside_s"]), s(res["freq_time_outside_s"]), p95(res["freq_time_outside_s"])),
    }
    print(agg)

#     # Build LaTeX (single-system row). Duplicate/merge for 9-bus & 39-bus later if you want side-by-side.
#     def fmt(x): 
#         if x != x: return r"--"
#         return f"{x:.2f}"
#     lines = []
#     for metric, (mean, std, p95v) in agg.items():
#         lines.append(f"{metric} & {fmt(mean)} $\\pm$ {fmt(std)} (p95 {fmt(p95v)}) \\\\")
#     table = rf"""
# \begin{table}[t]
# \centering
# \caption{{Policy quality for {system_label} (mean $\pm$ std across episodes; $p95$ in parentheses).}}
# \label{{tab:policy_quality_{system_label.replace(' ','_').lower()}}}
# \begin{tabular}{lc}
# \hline
# \textbf{{Metric}} & \textbf{{{system_label}}} \\
# \hline
# {chr(10).join(lines)}
# \hline
# \end{tabular}
# \end{table}
# """.strip()

#     Path(out_tex).write_text(table)
    return res, agg#, table

def plot_policy_quality_timeseries(res: pd.DataFrame, outdir: str | Path, title_prefix: str = ""):
    """
    Save three plots vs episode:
      - load_shed_mwh
      - freq_nadir_hz
      - freq_time_outside_s
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Ensure sorted by episode for a clean x-axis
    if "episode" not in res.columns:
        raise ValueError("Expected 'res' to have an 'episode' column.")
    # breakpoint()
    res_sorted = res.sort_values("episode")[["load_shed_mwh", "freq_nadir_hz", "freq_time_outside_s"]].rolling(10, min_periods=1).mean()

    plots = [
        ("load_shed_mwh",       "Load shed (MWh / episode)",            "policy_load_shed_mwh.png"),
        ("freq_nadir_hz",       "Frequency nadir (Hz)",                 "policy_freq_nadir_hz.png"),
        ("freq_time_outside_s", "Time outside [59.5, 60.5] Hz (s)",     "policy_freq_time_outside_s.png"),
    ]

    for col, ylabel, fname in plots:
        if col not in res_sorted.columns:
            print(f"[warn] Column '{col}' not found in res; skipping plot.")
            continue
        plt.figure()
        plt.plot(res_sorted.index.values.tolist(), res_sorted[col].values)
        plt.xlabel("Episode")
        plt.ylabel(ylabel)
        if title_prefix:
            plt.title(f"{title_prefix}: {ylabel} vs episode")
        else:
            plt.title(f"{ylabel} vs episode")
        plt.tight_layout()
        plt.savefig(outdir / fname, dpi=200)
        plt.close()

# ---------- System Efficiency (simulate) ----------
def aggregate_episode_stats(ep_csv: Path) -> Optional[pd.DataFrame]:
    df = read_csv(ep_csv)
    if df is None or df.empty:
        return None
    return df

def mean_std(series: pd.Series) -> Tuple[float, float]:
    if series is None or len(series) == 0:
        return float("nan"), 0.0
    m = float(series.mean())
    s = float(series.std(ddof=1)) if len(series) > 1 else 0.0
    return m, s

def summarize_system(df: pd.DataFrame, label: str) -> dict:
    out = {"label": label}
    # T_wall ≈ n_steps * avg_step_latency_ms / 1000
    if {"n_steps", "latency_ms_avg"}.issubset(df.columns):
        twall = df["n_steps"] * df["latency_ms_avg"] / 1000.0
        out["twall_s_mean"], out["twall_s_std"] = mean_std(twall)
    else:
        out["twall_s_mean"] = out["twall_s_std"] = float("nan")
    for col in ["latency_ms_avg", "latency_ms_p95", "peak_rss_mib", "delta_rss_mib"]:
        if col in df.columns:
            m, s = mean_std(df[col])
        else:
            m, s = float("nan"), 0.0
        out[f"{col}_mean"], out[f"{col}_std"] = m, s
    return out

def latex_table_system_efficiency_inverted(rows: List[dict]) -> str:
    systems = [r["label"] for r in rows]
    def fmt(x):
        try:
            if x != x or x is None:
                return "--"
            return f"{x:.2f}"
        except:
            return "--"
    metrics = [
        ("Wall-clock runtime $T_{\\text{wall}}$ (s)", "twall_s_mean", "twall_s_std"),
        ("Avg step latency $\\ell_{\\text{env}}$ (ms)", "latency_ms_avg_mean", "latency_ms_avg_std"),
        ("p95 step latency (ms)", "latency_ms_p95_mean", "latency_ms_p95_std"),
        ("Peak RSS (MiB)", "peak_rss_mib_mean", "peak_rss_mib_std"),
        ("$\\Delta$RSS (MiB)", "delta_rss_mib_mean", "delta_rss_mib_std"),
    ]
    header = "\\begin{table}[t]\n\\renewcommand{\\arraystretch}{1.3}\n\\centering\n"
    header += "\\caption{System efficiency: wall-clock runtime, step latency, and memory (mean $\\pm$ std across episodes).}\n"
    header += "\\label{tab:bus_perf}\n"
    header += "\\begin{tabular}{l" + "c"*len(systems) + "}\n\\hline\n"
    header += "\\textbf{Metric} & " + " & ".join(systems) + " \\\\\n\\hline\n"
    body = ""
    for title, mean_key, std_key in metrics:
        row_cells = []
        for r in rows:
            m, s = fmt(r.get(mean_key)), fmt(r.get(std_key, 0.0))
            row_cells.append(f"{m} $\\pm$ {s}")
        body += f"{title} & " + " & ".join(row_cells) + " \\\\\n"
    footer = "\\hline\n\\end{tabular}\n\\end{table}\n"
    return header + body + footer

def plot_step_latency_hist(sim_dir: Path, out_png: Path, title: str) -> bool:
    files = sorted(sim_dir.glob("steps_ep*.csv"))
    if not files:
        print(f"[warn] No step traces in {sim_dir}")
        return False
    lat_all = []
    for f in files:
        df = read_csv(f)
        if df is None or df.empty or "latency_ms" not in df.columns:
            continue
        lat_all.append(df["latency_ms"].values)
    if not lat_all:
        print(f"[warn] No latency_ms column found in step traces.")
        return False
    lat = np.concatenate(lat_all)
    plt.figure()
    plt.hist(lat, bins=50)
    plt.xlabel("GridPACK Step Latency (ms)")
    plt.ylabel("Count")
    # plt.title(title)
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()
    return True

def plot_cumulative_reward(train_sim_dir: Path, out_png: Path, title: str) -> bool:
    csvp = train_sim_dir / "steps.csv"
    df = read_csv(csvp)
    if df is None or df.empty:
        print(f"[warn] No reward CSV at {csvp}")
        return False
    for col in ["episode", "t", "reward"]:
        if col not in df.columns:
            print(f"[warn] Missing column '{col}' in {csvp}")
            return False
    
    cum_per_ep = df.groupby("episode")["reward"].sum()
    rolling_mean = cum_per_ep.rolling(10, min_periods=1).mean()
    
    # --- BEGIN: Add smooth band ---
    import numpy as np
    from scipy.interpolate import interp1d

    rolling_std = cum_per_ep.rolling(10, min_periods=1).std()
    # Smooth the std slightly to reduce jaggedness
    rolling_std_smooth = rolling_std.rolling(5, min_periods=1).mean()
    upper_bound = rolling_mean + rolling_std_smooth
    lower_bound = rolling_mean - rolling_std_smooth

    # Interpolate bounds for smooth fill
    orig_x = rolling_mean.index.values
    dense_x = np.linspace(orig_x.min(), orig_x.max(), num=200)

    def interp(x, y, dense_x):
        valid = ~np.isnan(y)
        if not np.any(valid):
            return np.full_like(dense_x, np.nan)
        f = interp1d(x[valid], y[valid], kind='linear', fill_value="extrapolate")
        return f(dense_x)

    interp_upper = interp(orig_x, upper_bound.values, dense_x)
    interp_lower = interp(orig_x, lower_bound.values, dense_x)
    # --- END: Add smooth band ---

    plt.figure()
    plt.plot(cum_per_ep.index.values, rolling_mean.values)  # original line unchanged

    # Fill smooth band
    plt.fill_between(dense_x, interp_lower, interp_upper, color='blue', alpha=0.3)

    plt.xlabel("Episode")
    plt.ylabel("Cumulative reward (10-ep MA)")
    # plt.title(title)
    if "input_39bus_IBR_fault" in str(train_sim_dir):
        plt.ylim(-2500, 0)
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()
    return True


# ---------- Training plots (voltage / frequency / reward) ----------
def find_columns(df: pd.DataFrame, prefixes: List[str]) -> List[str]:
    cols = []
    for p in prefixes:
        cols.extend([c for c in df.columns if c.startswith(p)])
    return sorted(set(cols))

def plot_train_voltage(df: pd.DataFrame, out_png: Path, title: str) -> bool:
    # Try common voltage column prefixes
    vcols = [c for c in df.columns if re.match("bus_.*_vn_pu", c)]
    if not vcols:
        print("[warn] No voltage columns found (e.g., bus_v_*)")
        return False
    # Pick last episode
    ep = int(df["episode"].max()) if "episode" in df.columns else None
    dff = df[df["episode"] == ep] if ep is not None else df
    # Mean across buses per time step
    # grp = dff.groupby("tick")[vcols].mean()
    plt.figure()
    plt.plot(dff["tick"], dff[vcols])
    plt.xlabel("t (Grid2Op step)")
    plt.ylabel("Bus voltage (p.u.)")
    # plt.title(title)
    # plt.xlim([3.95, 4.4])
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()
    return True

def plot_train_frequency(df: pd.DataFrame, out_png: Path, title: str) -> bool:
    fcols = [c for c in df.columns if re.match("bus_.*_frequency", c)]
    if not fcols:
        print("[warn] No frequency columns found (e.g., bus_f_*)")
        return False
    ep = int(df["episode"].min()) if "episode" in df.columns else None
    dff = df[df["episode"] == ep] if ep is not None else df
    plt.figure()
    plt.plot(dff["tick"], dff[fcols])
    plt.xlabel("t (Grid2Op step)")
    plt.ylabel("Bus frequency (Hz)")
    plt.xlim([3.9, 4.4])
    # plt.title(title)
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()
    return True

# ---------- Main ----------
def main():
    ap = argparse.ArgumentParser(description="Generate results (simulate + train) from a run folder.")
    ap.add_argument("--rundir", type=str, help="Path to run folder containing simulate/, train/, train_sim/, deploy/")
    ap.add_argument("--outdir", type=str, required=True, help="Directory to write figures and tables")
    args = ap.parse_args()

    run_dir = Path(args.rundir)
    suffix = os.path.basename(args.rundir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # ----- SIMULATE -----
    simulate_dir = run_dir / "simulate"
    train_dir = run_dir / "train"
    train_sim_dir = run_dir / "train_sim"
    deploy_dir = run_dir / "deploy"
    
    # ----- DEPLOY (inference) -----
    if deploy_dir.exists():
        ep_csv = deploy_dir / "system_efficiency_episodes.csv"
        ep_df = read_csv(ep_csv)
        if ep_df is not None and not ep_df.empty:
            # Summarize (supports multiple systems if you aggregate later)
            rows = []
            if "system" in ep_df.columns:
                for label, dfg in ep_df.groupby("system"):
                    rows.append(summarize_deploy_inference(dfg, str(label)))
            else:
                rows.append(summarize_deploy_inference(ep_df, suffix))

            tex = latex_table_deploy_inference(rows)
            (outdir / f"deploy_inference_{suffix}.tex").write_text(tex)

            # Per-episode plots
            plot_deploy_lines(ep_df, outdir, suffix)

            # Optional histograms if per-step columns exist in steps_ep*.csv
            plot_deploy_step_hists_if_available(deploy_dir, outdir, suffix)
        else:
            print(f"[info] deploy: missing or empty {ep_csv}")
    else:
        print("[info] deploy/ missing; skipping inference outputs.")

    if simulate_dir.exists():
        ep_csv = simulate_dir / "system_efficiency_episodes.csv"
        df = aggregate_episode_stats(ep_csv)
        if df is not None:
            rows = []
            if "system" in df.columns:
                for label, dfi in df.groupby("system"):
                    rows.append(summarize_system(dfi, str(label)))
            else:
                rows.append(summarize_system(df, "System"))
            tex = latex_table_system_efficiency_inverted(rows)
            (outdir / f"system_efficiency_{suffix}.tex").write_text(tex)
        else:
            print("[info] simulate: no episode stats CSV found.")

        plot_step_latency_hist(simulate_dir, outdir / f"step_latency_{suffix}.png",
                               "Step latency distribution (simulate)")

        # Reward plot (from train_sim/sim_all.csv), as requested to include under simulate outputs
        if train_sim_dir.exists():
            plot_cumulative_reward(train_dir, outdir / f"reward_{suffix}.png",
                                                  "Cumulative reward (10-ep MA)")
        else:
            print("[info] train_sim/ missing; skipping reward plot for simulate.")

    else:
        print("[info] simulate/ missing; skipping simulate outputs.")

    snappyfile = train_sim_dir / "sim_all.pq.snappy"
    csvfile = simulate_dir / "sim_all.csv"
    if os.path.exists(snappyfile):
        df = pd.read_parquet(snappyfile)
        if df is None or df.empty:
            print(f"[warn] No sim_all.pq.snappy for voltage at {snappyfile}")
            return False
    elif os.path.exists(csvfile):
        df = read_csv(csvfile)
    else:
        df = None
        print(f"[info] No sim_all.pq.snappy at {snappyfile}; skipping train voltage/reward/frequency plots.")

    if df is not None and not df.empty:
        # ----- TRAIN -----
        if train_sim_dir.exists():
            plot_train_voltage(df, outdir / f"voltage_{suffix}.png",
                            "Training: mean bus voltage (Episode 1)")
            plot_train_frequency(df, outdir / f"frequency_{suffix}.png",
                                "Training: mean bus frequency (Episode 1)")
            res, agg = summarize_policy_quality(df)
            plot_policy_quality_timeseries(res, outdir, title_prefix=suffix)

            evaluate_actions(df, outdir, no_op_idx=0, effect_window=3, pre_window=3)
        else:
            print("[info] train_sim/ missing; skipping train voltage/reward/frequency plots.")
    
    print(f"[done] Wrote outputs to: {outdir}")

if __name__ == "__main__":
    main()
