"""NA needed to reach a target number of resolvable synapses, vs density.

For each labeling density (0.01-0.30, 15 pts) and NA (0.30-1.00), the resolvable
synapse count (purity > threshold) is summed across 10 random FOVs (run in
parallel). Then for each target count (1000, 2000, 4000) we find the smallest NA
that reaches it, and plot NA-needed vs density - one colored curve per target.

    python scripts/fig_na_for_counts.py
"""
import os
from dataclasses import replace
from multiprocessing import Pool

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from synsim import DEFAULTS, build_scene, evaluate

DENSITIES = np.linspace(0.01, 0.30, 15)
NA_GRID = np.round(np.arange(0.30, 1.001, 0.05), 3)
REGIONS = list(range(1, 11))
TARGETS = [(1000, "tab:blue"), (2000, "tab:orange"), (4000, "tab:green")]


def _work(region):
    """One FOV -> resolvable-count matrix [density, NA]."""
    scene = build_scene(rng=np.random.default_rng(0),
                        params=replace(DEFAULTS, region=float(region)))
    u = np.random.default_rng(100 + region).random(len(scene.axons))
    counts = np.zeros((len(DENSITIES), len(NA_GRID)))
    for di, dens in enumerate(DENSITIES):
        mask = u < dens
        for ni, na in enumerate(NA_GRID):
            counts[di, ni] = evaluate(scene, mask, replace(DEFAULTS, na=float(na)))["resolvable"]
    return counts


def first_crossing(na_grid, vals, target):
    """Smallest NA at which vals reaches target (first up-crossing, lin. interp)."""
    vals = np.asarray(vals)
    above = vals >= target
    if not above.any():
        return np.nan
    i = int(np.argmax(above))
    if i == 0:
        return na_grid[0]
    x0, x1, y0, y1 = na_grid[i - 1], na_grid[i], vals[i - 1], vals[i]
    return x1 if y1 == y0 else x0 + (target - y0) * (x1 - x0) / (y1 - y0)


def main():
    nproc = min(len(REGIONS), os.cpu_count() or 1)
    print(f"running {len(REGIONS)} FOVs on {nproc} cores...", flush=True)
    with Pool(nproc) as pool:
        mats = pool.map(_work, REGIONS)
    total = np.sum(mats, axis=0)            # [density, NA] pooled resolvable count
    imax = np.unravel_index(int(np.argmax(total)), total.shape)
    print(f"peak pooled resolvable = {total.max():.0f} at density "
          f"{DENSITIES[imax[0]]:.3f}, NA {NA_GRID[imax[1]]:.2f}")

    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    print("\n density   " + "   ".join(f"NA@{t}" for t, _ in TARGETS))
    unreachable = []
    for target, color in TARGETS:
        na_need = np.array([first_crossing(NA_GRID, total[di], target)
                            for di in range(len(DENSITIES))])
        ax.plot(DENSITIES, na_need, "o-", color=color, label=f"{target} synapses")
        if np.isnan(na_need).all():
            unreachable.append((target, color))
    for di, d in enumerate(DENSITIES):
        vals = [first_crossing(NA_GRID, total[di], t) for t, _ in TARGETS]
        print(f" {d:5.3f}    " + "   ".join(f"{v:.3f}" if not np.isnan(v) else "  -  " for v in vals))

    ax.axhline(0.6, color="0.5", ls=":", lw=1)
    ax.annotate("mesoscope NA 0.6", (DENSITIES[-1], 0.6), fontsize=8, va="bottom", ha="right", color="0.4")
    if unreachable:
        labels = ", ".join(str(t) for t, _ in unreachable)
        ax.text(0.5, 0.4, f"{labels} synapses: not reachable\n(peak yield {total.max():.0f} "
                f"over 10 FOVs at density {DENSITIES[imax[0]]:.2f}, NA {NA_GRID[imax[1]]:.2f})",
                transform=ax.transAxes, ha="center", fontsize=9, color="tab:green",
                bbox=dict(boxstyle="round", fc="white", ec="tab:green"))
    ax.set_xlabel("labeling density (fraction of axons)")
    ax.set_ylabel("NA needed")
    ax.set_ylim(0.3, 1.05)
    ax.set_title("NA needed to reach N resolvable synapses (total over 10 FOVs ~ 30x30x19 um each)")
    ax.legend(title="target")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig("figures/na_for_counts.png", dpi=120)
    print("-> figures/na_for_counts.png")


if __name__ == "__main__":
    main()
