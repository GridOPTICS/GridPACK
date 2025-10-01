#!/usr/bin/env python3
from pathlib import Path
import numpy as np
import pandas as pd

# --------- helpers to autodetect columns ---------
def _cols(df, suffix):
    return [c for c in df.columns if c.endswith(suffix)]

def bus_voltage_cols(df):     return _cols(df, "_vn_pu")
def bus_frequency_cols(df):   return _cols(df, "_frequency")
def load_p_cols(df):          return [c for c in df.columns if c.startswith("load_") and c.endswith("_p_mw")]

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
def summarize_policy_quality(sim_all_pq: str, out_tex: str, system_label: str):
    df = pd.read_parquet(sim_all_pq)
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
    breakpoint()

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

if __name__ == "__main__":
    # --------- example usage ---------
    res9, agg9 = summarize_policy_quality(
        "runs/dqn_grid2op/input_9bus_fault/train_sim/sim_all.pq.snappy",
        "results/table_policy_quality_9bus.tex",
        "IEEE 9-bus"
    )
    res9.to_csv("./runs/policy_quality_9bus.csv", index=False)
    # res39, agg39, tex39 = summarize_policy_quality("runs/39bus/train_sim/sim_all.csv",
    #                                               "results/table_policy_quality_39bus.tex",
    #                                               "IEEE 39-bus")
