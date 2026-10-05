"""Shared plotting for the power-curve figures."""

SERIES, INK, INK2, MUTED, GRID, AXIS, SURFACE = (
    "#2a78d6", "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7", "#fcfcfb")
CATEGORICAL = ["#2a78d6", "#eb6834"]   # fixed order: slot 1 blue, slot 2 orange


def plot_curves(x, xlabel, panels, title, path):
    """Side-by-side line charts sharing an x-axis.

    panels: list of dicts with keys
      y      -- values to plot against x
      series -- instead of y: list of (label, values) pairs, drawn in the
                CATEGORICAL colors in order, with a legend
      legend_loc -- optional matplotlib legend location for series panels
      title  -- panel title
      ylim   -- (lo, hi), or None for 0 to just above the largest value
      ylabel -- optional y-axis label
      ref    -- optional (value, label) for a horizontal reference line
      mark_x -- optional (x value, label) for a vertical reference line
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, len(panels), figsize=(5 * len(panels), 4.2),
                             facecolor=SURFACE, squeeze=False)
    for ax, p in zip(axes.ravel(), panels):
        ax.set_facecolor(SURFACE)
        series = p.get("series") or [(None, p["y"])]
        for (label, y), color in zip(series, CATEGORICAL if len(series) > 1 else [SERIES]):
            ax.plot(x, y, color=color, lw=2, solid_joinstyle="round",
                    solid_capstyle="round", marker="o", ms=8,
                    markeredgecolor=SURFACE, markeredgewidth=2, zorder=3, label=label)
        if len(series) > 1:
            ax.legend(loc=p.get("legend_loc", "best"), frameon=False, labelcolor=INK2,
                      fontsize=9)
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
        top = max(max(y) for _, y in series)
        ax.set_ylim(*(p.get("ylim") or (0, top * 1.08)))
        if p.get("ylabel"):
            ax.set_ylabel(p["ylabel"], color=INK2)
        if p.get("ref"):
            value, label = p["ref"]
            lo, hi = ax.get_ylim()
            ax.axhline(value, color=MUTED, lw=1, zorder=2)
            ax.text(x[0], value + 0.02 * (hi - lo), label, color=INK2, fontsize=9)
        if p.get("mark_x"):
            value, label = p["mark_x"]
            lo, hi = ax.get_ylim()
            ax.axvline(value, color=MUTED, lw=1, ls=(0, (4, 3)), zorder=2)
            ax.text(value, lo + 0.02 * (hi - lo), f"{label} ", color=INK2, fontsize=9,
                    ha="right", va="bottom")
    fig.suptitle(title, x=0.01, ha="left", color=INK2, fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=150, facecolor=SURFACE)
    plt.close(fig)
