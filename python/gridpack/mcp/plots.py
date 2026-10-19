# -------------------------------------------------------------
# file: gridpack/mcp/plots.py
# -------------------------------------------------------------
# Plots of job results (needs matplotlib: gridpack[plot]).
# -------------------------------------------------------------

"""Render job results as a figure: a PNG for the chat, or a PDF, SVG or
PNG file in a style fit for a paper.  Imports nothing from gridpack so it
can be exercised without the extension; matplotlib is loaded lazily.
"""

import fnmatch
import io

__all__ = ["FORMATS", "STYLES", "WIDTHS", "build", "render", "save"]

# Categorical palette: fixed order, never cycled past eight series.
PALETTE = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300",
           "#6250d6", "#e34948"]
# Dash cycle so a paper figure still reads in grayscale.
DASHES = ["-", "--", "-.", ":", (0, (5, 1, 1, 1)), (0, (3, 1, 1, 1, 1, 1)),
          (0, (1, 1)), (0, (6, 2))]
_INK, _MUTED, _GRID, _BAND = "#15181c", "#5b6470", "#e3e6ea", "#15181c"
_UNITS = {"V": "p.u.", "Pg": "p.u.", "Qg": "p.u.", "angle": "rad", "speed": "p.u.",
          "Efd": "p.u.", "Pm": "p.u.", "PowerAngle": "deg or rad (model)",
          "rspeed": "p.u.", "rangle": "rad", "P": "p.u.", "Q": "p.u.",
          "vmag": "p.u.", "vangle": "rad", "freq": "Hz", "online": "fraction"}
_LABELS = {"speed": "rotor speed", "rspeed": "rotor speed", "angle": "rotor angle",
           "rangle": "rotor angle", "V": "terminal voltage", "vmag": "bus voltage",
           "vangle": "bus angle", "Pg": "electrical power", "P": "electrical power",
           "Qg": "reactive power", "Q": "reactive power", "Efd": "field voltage",
           "Pm": "mechanical power", "freq": "bus frequency",
           "PowerAngle": "power angle"}
_INSTALL = "plot_job needs matplotlib: pip install 'gridpack[plot]'"

FORMATS = ("png", "pdf", "svg")
STYLES = ("screen", "paper")
# Column widths in inches: IEEE single and double column.
WIDTHS = {"single": 3.5, "double": 7.16}


def _mpl(style, font):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        raise ValueError(_INSTALL)
    rc = {"text.color": _INK, "axes.labelcolor": _MUTED, "xtick.color": _MUTED,
          "ytick.color": _MUTED, "axes.edgecolor": _GRID, "axes.titlecolor": _INK,
          "font.family": "serif" if font == "serif" else "sans-serif",
          "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none"}
    if style == "paper":
        rc.update({"font.size": 8, "axes.labelsize": 8, "xtick.labelsize": 7,
                   "ytick.labelsize": 7, "legend.fontsize": 7,
                   "axes.labelcolor": _INK, "xtick.color": _INK,
                   "ytick.color": _INK, "axes.edgecolor": _INK,
                   "lines.linewidth": 1.0})
    else:
        rc.update({"font.size": 10, "lines.linewidth": 2.0})
    plt.rcParams.update(rc)
    return plt


def render(result, app, kind=None, channels=None, title=None, events=()):
    """PNG bytes for the chat (screen style)."""
    fig, plt = build([result], [app], kind, channels, title, [events])
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=150, facecolor="white")
    plt.close(fig)
    return buf.getvalue()


def save(path, results, apps, kind=None, channels=None, title=None,
         events_list=None, labels=None, style="paper", fmt="pdf",
         width=None, height=None, font="serif"):
    """Write the figure to ``path``; return (width, height) in inches."""
    if fmt not in FORMATS:
        raise ValueError("format is one of %s, not %r" % (", ".join(FORMATS), fmt))
    fig, plt = build(results, apps, kind, channels, title, events_list, labels,
                     style, width, height, font)
    size = tuple(float(x) for x in fig.get_size_inches())
    fig.savefig(str(path), format=fmt, dpi=300 if fmt == "png" else None,
                facecolor="white", bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
    return size


def build(results, apps, kind=None, channels=None, title=None, events_list=None,
          labels=None, style="screen", width=None, height=None, font="sans"):
    """The matplotlib figure for one or more jobs of the same kind.

    Dynamic runs: ``kind`` "series" over ``channels`` (names or globs;
    default every speed), one panel per quantity.  Several jobs overlay:
    one color per job, one dash per channel.  Power flow / state
    estimation: "voltage" or "angle" profiles over the buses.  Paper style
    drops the title, uses print sizes and dashes every series.
    """
    if style not in STYLES:
        raise ValueError("style is one of %s, not %r" % (", ".join(STYLES), style))
    if not results or len(results) != len(apps):
        raise ValueError("give one result per job")
    if len(set("dyn" if a in ("dsf", "hadrec") else a for a in apps)) > 1:
        raise ValueError("all jobs must be of the same kind; got %s" % ", ".join(apps))
    plt = _mpl(style, font)
    w = _width(width)
    labels = list(labels or ["job %d" % (i + 1) for i in range(len(results))])
    if len(labels) != len(results):
        raise ValueError("give one label per job")
    if apps[0] in ("dsf", "hadrec"):
        fig = _series(plt, results, apps, kind or "series", channels, title,
                      events_list or [()] * len(results), labels, style, w, height)
    elif apps[0] in ("powerflow", "se"):
        fig = _profile(plt, results, kind or "voltage", title, labels, style, w, height)
    else:
        raise ValueError("no plot for a %s job; use job_report" % apps[0])
    return fig, plt


def _width(width):
    if width is None:
        return None
    if width in WIDTHS:
        return WIDTHS[width]
    if isinstance(width, bool) or not isinstance(width, (int, float)) or width <= 0:
        raise ValueError("width is %s or inches" % " or ".join(WIDTHS))
    return float(width)


def _quantity(col):
    return col.rsplit("_", 1)[-1]


def _pick(cols, channels):
    picked = []
    for p in channels or ["*_speed", "*_rspeed"]:
        hits = [c for c in cols[1:] if c == p or fnmatch.fnmatchcase(c, p)]
        if not hits and channels:
            raise ValueError("no channel matches %r; channels are %s"
                             % (p, ", ".join(cols[1:])))
        picked += [h for h in hits if h not in picked]
    if not picked:
        raise ValueError("no speed channel; give channels, from %s" % ", ".join(cols[1:]))
    return picked


def _series(plt, results, apps, kind, channels, title, events_list, labels,
            style, width, height):
    if kind != "series":
        raise ValueError("a %s job plots kind \"series\", not %r" % (apps[0], kind))
    tables = []
    for r, a in zip(results, apps):
        data = r.get("watch" if a == "dsf" else "observations")
        if not data or not data.get("rows"):
            raise ValueError("a job recorded no time series")
        tables.append((data["columns"], data["rows"]))
    cols = tables[0][0]
    picked = _pick(cols, channels)
    for c, _ in tables[1:]:
        missing = [p for p in picked if p not in c]
        if missing:
            raise ValueError("not every job has channel %s" % ", ".join(missing))
    groups = {}
    for c in picked:
        groups.setdefault(_quantity(c), []).append(c)
    many = len(results) > 1
    for q, cs in groups.items():
        if (len(cs) if not many else len(results)) > len(PALETTE):
            raise ValueError("too many series for one panel (%d); pick at most %d "
                             "channels or jobs" % (max(len(cs), len(results)),
                                                   len(PALETTE)))
    n = len(groups)
    w = width or 9.0
    h = height or (n * (0.62 if style == "paper" else 0.29) * w + (0.3 if style == "paper" else 0.8))
    fig, axes = plt.subplots(n, 1, figsize=(w, h), sharex=True, squeeze=False)
    for ax, (q, cs) in zip(axes[:, 0], groups.items()):
        for j, ((tcols, rows), lab) in enumerate(zip(tables, labels)):
            t = [row[0] for row in rows]
            for k, c in enumerate(cs):
                i = tcols.index(c)
                color = PALETTE[j if many else k]
                # Grayscale safety: dash by channel, or by job when the
                # jobs share one channel; screen style keeps solid lines
                # for a single job.
                if many:
                    dash = DASHES[k if len(cs) > 1 else j]
                else:
                    dash = DASHES[k if style == "paper" else 0]
                name = "%s %s" % (lab, c) if many and len(cs) > 1 else (lab if many else c)
                ax.plot(t, [row[i] for row in rows], color=color, linestyle=dash,
                        solid_joinstyle="round", label=name)
        _events(ax, events_list[0])
        ax.set_ylabel("%s (%s)" % (_LABELS.get(q, q), _UNITS.get(q, "")))
        _axes(ax, style)
        ax.set_xlim(min(row[0] for row in tables[0][1]),
                    max(row[0] for row in tables[0][1]))
        if many or len(cs) > 1:
            ax.legend(frameon=False, loc="best", ncol=1 if style == "paper" else 2)
    axes[-1, 0].set_xlabel("time (s)")
    if style != "paper":
        axes[0, 0].set_title(title or ("%s run%s, %d channel%s" % (
            "Dynamic simulation" if apps[0] == "dsf" else "HADREC",
            "s" if many else "", len(picked), "s" if len(picked) != 1 else "")),
            loc="left", fontweight="bold")
    elif title:
        axes[0, 0].set_title(title, loc="left")
    fig.tight_layout()
    return fig


def _events(ax, events):
    for e in events or ():
        if "start" in e:
            ax.axvspan(e["start"], e["end"], color=_BAND, alpha=0.07, lw=0)
        else:
            ax.axvline(e["time"], color=_MUTED, lw=0.8, ls=(0, (3, 3)))


def _axes(ax, style):
    ax.grid(axis="y", color=_GRID, lw=0.6 if style == "paper" else 1)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.tick_params(length=2 if style == "paper" else 0)


def _profile(plt, results, kind, title, labels, style, width, height):
    key = {"voltage": "voltage", "angle": "angle"}.get(kind)
    if key is None:
        raise ValueError("kind is \"voltage\" or \"angle\", not %r" % kind)
    tables = []
    for r in results:
        buses = r.get("buses", [])
        if not buses:
            raise ValueError("a job has no bus table")
        tables.append(sorted(buses, key=lambda b: b["busId"]))
    ids = [b["busId"] for b in tables[0]]
    w = width or 9.0
    h = height or ((0.55 if style == "paper" else 0.4) * w)
    fig, ax = plt.subplots(figsize=(w, h))
    many = len(results) > 1
    for j, (buses, lab) in enumerate(zip(tables, labels)):
        y = {b["busId"]: b[key] for b in buses}
        ax.plot(range(len(ids)), [y.get(i, float("nan")) for i in ids],
                color=PALETTE[j], linestyle=DASHES[j if style == "paper" else 0],
                marker="o" if len(ids) <= 60 else None,
                ms=(2.5 if style == "paper" else 4), label=lab)
    if key == "voltage":
        for lim in (0.9, 1.1):
            ax.axhline(lim, color=_MUTED, lw=0.8, ls=(0, (3, 3)))
        ax.set_ylabel("voltage (p.u.)")
    else:
        ax.set_ylabel("angle (deg)")
    if len(ids) <= 30:
        ax.set_xticks(range(len(ids)))
        ax.set_xticklabels([str(i) for i in ids], rotation=90 if len(ids) > 14 else 0)
        ax.set_xlabel("bus")
    else:
        ax.set_xlabel("bus (in number order, %d buses)" % len(ids))
    _axes(ax, style)
    if many:
        ax.legend(frameon=False, loc="best")
    if style != "paper":
        ax.set_title(title or "Bus %s profile" % key, loc="left", fontweight="bold")
    elif title:
        ax.set_title(title, loc="left")
    fig.tight_layout()
    return fig
