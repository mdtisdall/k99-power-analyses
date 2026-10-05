"""Shared plotting for the power-curve figures in equivalence.py and group_effect.py."""

SERIES, INK, INK2, MUTED, GRID, AXIS, SURFACE = (
    "#2a78d6", "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7", "#fcfcfb")


def plot_curves(x, xlabel, panels, title, path):
    """Side-by-side single-series line charts sharing an x-axis.

    panels: list of dicts with keys
      y      -- values to plot against x
      title  -- panel title
      ylim   -- (lo, hi), or None for 0 to just above the largest value
      ref    -- optional (value, label) for a horizontal reference line
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, len(panels), figsize=(5 * len(panels), 4.2),
                             facecolor=SURFACE)
    for ax, p in zip(axes, panels):
        ax.set_facecolor(SURFACE)
        ax.plot(x, p["y"], color=SERIES, lw=2, solid_joinstyle="round",
                solid_capstyle="round", marker="o", ms=8,
                markeredgecolor=SURFACE, markeredgewidth=2, zorder=3)
        ax.set_title(p["title"], loc="left", color=INK, fontsize=11)
        ax.set_xlabel(xlabel, color=INK2)
        # Label every x value unless they are too close together to read.
        gaps = [b - a for a, b in zip(x, x[1:])]
        if not gaps or min(gaps) >= 0.06 * (x[-1] - x[0]):
            ax.set_xticks(x)
        ax.grid(True, axis="y", color=GRID, lw=1)
        ax.set_axisbelow(True)
        for side in ("top", "right", "left"):
            ax.spines[side].set_visible(False)
        ax.spines["bottom"].set_color(AXIS)
        ax.tick_params(colors=MUTED, length=0)
        ax.set_ylim(*(p.get("ylim") or (0, max(p["y"]) * 1.08)))
        if p.get("ref"):
            value, label = p["ref"]
            lo, hi = ax.get_ylim()
            ax.axhline(value, color=MUTED, lw=1, zorder=2)
            ax.text(x[0], value + 0.02 * (hi - lo), label, color=INK2, fontsize=9)
    fig.suptitle(title, x=0.01, ha="left", color=INK2, fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=150, facecolor=SURFACE)
    plt.close(fig)
