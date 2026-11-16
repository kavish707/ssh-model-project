"""Shared figure style for Project A.

Categorical hues are used in a fixed order (never cycled) and were validated with
the dataviz palette checker (lightness band, chroma floor, CVD separation,
normal-vision floor all PASS in light mode on surface #fcfcfb).  Because three of
the hues sit below 3:1 contrast against the surface, every figure carries either a
legend with text labels or direct labels, and every numeric result is also written
to CSV -- the "table view" relief.  Line style and marker shape duplicate the
colour information, so identity is never carried by colour alone.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

SURFACE = "#fcfcfb"
INK = "#1c1c1a"
INK_SECONDARY = "#55554f"
INK_MUTED = "#8a8a82"
GRID = "#e3e3de"

# fixed categorical order -- do not shuffle, do not cycle
CAT = ["#0072B2", "#D55E00", "#009E73", "#E69F00", "#CC79A7", "#56B4E9"]
LINESTYLES = ["-", "--", "-.", (0, (3, 1, 1, 1)), ":", (0, (5, 1))]
MARKERS = ["o", "s", "^", "D", "v", "P"]

REFERENCE = INK_MUTED  # analytic / exact reference curves


def use_style() -> None:
    plt.rcParams.update(
        {
            "figure.facecolor": SURFACE,
            "axes.facecolor": SURFACE,
            "savefig.facecolor": SURFACE,
            "axes.edgecolor": INK_MUTED,
            "axes.labelcolor": INK,
            "axes.titlecolor": INK,
            "axes.linewidth": 0.8,
            "axes.grid": True,
            "axes.axisbelow": True,
            "grid.color": GRID,
            "grid.linewidth": 0.6,
            "xtick.color": INK_SECONDARY,
            "ytick.color": INK_SECONDARY,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "axes.labelsize": 10,
            "axes.titlesize": 10.5,
            "legend.fontsize": 8.5,
            "legend.frameon": False,
            "lines.linewidth": 1.8,
            "lines.markersize": 4.5,
            "font.size": 10,
            "figure.dpi": 120,
            "savefig.dpi": 200,
            "savefig.bbox": "tight",
            "mathtext.fontset": "dejavusans",
        }
    )


def style_axes(ax, title=None, xlabel=None, ylabel=None):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    if title:
        ax.set_title(title, loc="left", pad=8)
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    return ax
