# -------------------------------------------------------------
# file: gridpack/mcp/plots.py
# -------------------------------------------------------------
# PNG plots of job results (needs matplotlib: gridpack[plot]).
# -------------------------------------------------------------

"""Render a job's results to PNG bytes.  Imports nothing from gridpack so
it can be exercised without the extension; matplotlib is loaded lazily.
"""

import fnmatch
import io
import math

__all__ = ["render"]

# Categorical palette: fixed order, never cycled past eight series.
PALETTE = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300",
           "#6250d6", "#e34948"]
_INK, _MUTED, _GRID, _BAND = "#15181c", "#5b6470", "#e3e6ea", "#15181c"
_UNITS = {"V": "p.u.", "Pg": "p.u.", "Qg": "p.u.", "angle": "rad", "speed": "p.u.",
          "Efd": "p.u.", "Pm": "p.u.", "PowerAngle": "rad", "rspeed": "p.u.",
          "rangle": "rad", "P": "p.u.", "Q": "p.u.", "vmag": "p.u.",
          "vangle": "rad", "freq": "Hz", "online": "fraction"}
_INSTALL = "plot_job needs matplotlib: pip install 'gridpack[plot]'"


def _mpl():
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        raise ValueError(_INSTALL)
    plt.rcParams.update({"font.size": 10, "text.color": _INK,
                         "axes.labelcolor": _MUTED, "xtick.color": _MUTED,
                         "ytick.color": _MUTED, "axes.edgecolor": _GRID,
                         "axes.titlecolor": _INK})
    return plt


def render(result, app, kind=None, channels=None, title=None, events=()):
    """PNG bytes.  Dynamic runs: ``kind`` "series" (default) over
    ``channels`` (names or globs; default every speed), one panel per
    quantity.  Power flow / state estimation: "voltage" (default) or
    "angle" profiles over the buses."""
    plt = _mpl()
    if app in ("dsf", "hadrec"):
        fig = _series(plt, result, app, kind or "series", channels, title, events)
    elif app in ("powerflow", "se"):
        fig = _profile(plt, result, kind or "voltage", title)
    else:
        raise ValueError("no plot for a %s job; use job_report" % app)
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=150, facecolor="white")
    plt.close(fig)
    return buf.getvalue()


def _quantity(col):
    return col.rsplit("_", 1)[-1]


def _series(plt, result, app, kind, channels, title, events):
    if kind != "series":
        raise ValueError("a %s job plots kind \"series\", not %r" % (app, kind))
    data = result.get("watch" if app == "dsf" else "observations")
    if not data or not data.get("rows"):
        raise ValueError("the job recorded no time series")
    cols, rows = data["columns"], data["rows"]
    picked = []
    for p in channels or ["*_speed", "*_rspeed"]:
        hits = [c for c in cols[1:] if c == p or fnmatch.fnmatchcase(c, p)]
        if not hits and channels:
            raise ValueError("no channel matches %r; channels are %s"
                             % (p, ", ".join(cols[1:])))
        picked += [h for h in hits if h not in picked]
    if not picked:
        raise ValueError("no speed channel; give channels, from %s" % ", ".join(cols[1:]))
    groups = {}
    for c in picked:
        groups.setdefault(_quantity(c), []).append(c)
    for q, cs in groups.items():
        if len(cs) > len(PALETTE):
            raise ValueError("%d %s channels is too many for one plot; pick at "
                             "most %d" % (len(cs), q, len(PALETTE)))
    t = [r[0] for r in rows]
    fig, axes = plt.subplots(len(groups), 1, figsize=(9, 2.6 * len(groups) + 0.8),
                             sharex=True, squeeze=False)
    for ax, (q, cs) in zip(axes[:, 0], groups.items()):
        for k, c in enumerate(cs):
            i = cols.index(c)
            ax.plot(t, [r[i] for r in rows], color=PALETTE[k], lw=2,
                    solid_joinstyle="round", label=c)
        _events(ax, events)
        ax.set_ylabel("%s (%s)" % (q, _UNITS.get(q, "")))
        ax.grid(axis="y", color=_GRID, lw=1)
        ax.set_axisbelow(True)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        ax.tick_params(length=0)
        ax.set_xlim(t[0], t[-1])
        if len(cs) > 1:
            ax.legend(frameon=False, fontsize=9, loc="best")
    axes[-1, 0].set_xlabel("time (s)")
    axes[0, 0].set_title(title or ("%s run, %d channel%s" % (
        "Dynamic simulation" if app == "dsf" else "HADREC", len(picked),
        "s" if len(picked) != 1 else "")), loc="left", fontweight="bold")
    fig.tight_layout()
    return fig


def _events(ax, events):
    for e in events or ():
        if "start" in e:
            ax.axvspan(e["start"], e["end"], color=_BAND, alpha=0.07, lw=0)
        else:
            ax.axvline(e["time"], color=_MUTED, lw=1, ls=(0, (3, 3)))


def _profile(plt, result, kind, title):
    buses = result.get("buses", [])
    if not buses:
        raise ValueError("the job has no bus table")
    key = {"voltage": "voltage", "angle": "angle"}.get(kind)
    if key is None:
        raise ValueError("kind is \"voltage\" or \"angle\", not %r" % kind)
    buses = sorted(buses, key=lambda b: b["busId"])
    x = list(range(len(buses)))
    y = [b[key] for b in buses]
    fig, ax = plt.subplots(figsize=(9, 3.6))
    ax.plot(x, y, color=PALETTE[0], lw=2, marker="o", ms=4 if len(x) <= 60 else 0)
    if key == "voltage":
        for lim in (0.9, 1.1):
            ax.axhline(lim, color=_MUTED, lw=1, ls=(0, (3, 3)))
        ax.set_ylabel("voltage (p.u.)")
    else:
        ax.set_ylabel("angle (deg)")
    if len(x) <= 30:
        ax.set_xticks(x)
        ax.set_xticklabels([str(b["busId"]) for b in buses], rotation=90 if len(x) > 14 else 0)
        ax.set_xlabel("bus")
    else:
        ax.set_xlabel("bus (in number order, %d buses)" % len(x))
    ax.grid(axis="y", color=_GRID, lw=1)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.tick_params(length=0)
    ax.set_title(title or "Bus %s profile" % key, loc="left", fontweight="bold")
    fig.tight_layout()
    return fig
