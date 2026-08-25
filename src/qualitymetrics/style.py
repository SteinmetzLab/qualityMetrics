"""Lab plotting conventions, applied once so every figure obeys them.

Arial, no top/right spines, editable text in vector output, sentence-case
labels with units. Importing :func:`use_lab_style` and calling it is the only
thing a plotting module has to do; nothing here reaches into a figure after the
fact, so a caller who wants something different can simply not call it.
"""
from __future__ import annotations

import matplotlib as mpl
import numpy as np

#: Firing rate is spikes per second. "Hz" invites confusion with the
#: frequency-domain quantities that appear all over the other QC figures.
RATE_LABEL = "Firing rate (spikes/s)"
DEPTH_LABEL = "Depth (\u00b5m)"
TIME_LABEL = "Time (s)"
AMP_LABEL = "Amplitude (\u00b5V)"


def use_lab_style() -> None:
    """Set the rcParams the lab uses for every figure."""
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "axes.spines.top": False,
        "axes.spines.right": False,
        # Type 42 keeps text as text in PDF/PS, so a figure stays editable in
        # Illustrator instead of arriving as outlines.
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "figure.dpi": 110,
        "savefig.dpi": 200,
        "savefig.bbox": "tight",
        "axes.titlesize": 11,
        "axes.labelsize": 10,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "legend.fontsize": 9,
        "legend.frameon": False,
        "image.cmap": "viridis",
    })


def despine(ax) -> None:
    """Hide the top and right spines of an axes that was built by hand."""
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def size_legend(ax, sizes, labels, title, **kw):
    """A legend that explains what marker *area* means.

    The MATLAB original scaled marker size by spike count and never said so,
    which makes the most eye-catching channel in the figure the one piece of
    information a reader cannot decode. Any plot that varies marker size should
    call this.
    """
    from matplotlib.lines import Line2D
    handles = [Line2D([], [], marker="o", linestyle="none", color="0.35",
                      markersize=np.sqrt(s), label=str(lab))
               for s, lab in zip(sizes, labels, strict=False)]
    leg = ax.legend(handles=handles, title=title, labelspacing=1.1,
                    borderpad=0.8, handletextpad=1.0, **kw)
    leg.get_title().set_fontsize(9)
    return leg


def save(fig, path, close: bool = True) -> str:
    """Write a figure with a white ground, and report where it went."""
    fig.patch.set_facecolor("white")
    fig.savefig(path, facecolor="white")
    if close:
        import matplotlib.pyplot as plt
        plt.close(fig)
    return str(path)
