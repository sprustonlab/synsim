"""Three figures showing the advantage of higher NA (resolvable synapse yield).

Computes pooled resolvable count over 10 FOVs across density x NA (parallel), then:
  figures/yield_vs_density_byNA.png  - yield vs density, one curve per NA (family)
  figures/peak_yield_vs_na.png       - peak yield vs NA + optimal density (twin axis)
  figures/yield_heatmap.png          - resolvable count over density x NA, with contours

    python scripts/fig_na_advantage.py
"""
import os
from dataclasses import replace
from multiprocessing import Pool

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from synsim import DEFAULTS, build_scene, evaluate

DENSITIES = np.linspace(0.01, 0.60, 20)
NA_GRID = np.round(np.arange(0.30, 1.001, 0.05), 3)
REGIONS = list(range(1, 11))
NA_SHOW = [0.4, 0.5, 0.6, 0.7, 0.8, 0.9]


def _work(region):
    scene = build_scene(rng=np.random.default_rng(0),
                        params=replace(DEFAULTS, region=float(region)))
    u = np.random.default_rng(100 + region).random(len(scene.axons))
    counts = np.zeros((len(DENSITIES), len(NA_GRID)))
    for di, dens in enumerate(DENSITIES):
        mask = u < dens
        for ni, na in enumerate(NA_GRID):
            counts[di, ni] = evaluate(scene, mask, replace(DEFAULTS, na=float(na)))["resolvable"]
    return counts


def main():
    nproc = min(len(REGIONS), os.cpu_count() or 1)
    print(f"running {len(REGIONS)} FOVs on {nproc} cores...", flush=True)
    with Pool(nproc) as pool:
        mats = pool.map(_work, REGIONS)
    total = np.sum(mats, axis=0)            # [density, NA]

    # ---- Fig A: yield vs density, one curve per NA ----
    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    for na in NA_SHOW:
        ni = int(np.argmin(np.abs(NA_GRID - na)))
        col = plt.cm.viridis((na - 0.4) / 0.5)
        y = total[:, ni]
        ax.plot(DENSITIES, y, "o-", color=col, label=f"NA {na}")
        di = int(np.argmax(y))
        ax.scatter([DENSITIES[di]], [y[di]], s=90, facecolor="none", edgecolor=col, lw=2, zorder=5)
    ax.set_xlabel("labeling density (fraction of axons)")
    ax.set_ylabel("resolvable synapses (total over 10 FOVs)")
    ax.set_title("Higher NA raises the yield ceiling and shifts the optimal density right\n"
                 "(circles = peak per NA)")
    ax.legend(title="NA"); ax.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig("figures/yield_vs_density_byNA.png", dpi=120)

    # ---- Fig B: peak yield vs NA + optimal density ----
    peak = total.max(axis=0)
    opt_density = DENSITIES[np.argmax(total, axis=0)]
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(NA_GRID, peak, "o-", color="tab:blue")
    ax.set_xlabel("NA"); ax.set_ylabel("peak resolvable synapses (10 FOVs)", color="tab:blue")
    ax.tick_params(axis="y", labelcolor="tab:blue")
    ax.axvline(0.6, color="0.6", ls=":", lw=1); ax.annotate("mesoscope 0.6", (0.6, peak.min()),
                                                            fontsize=8, color="0.4", rotation=90, va="bottom")
    ax2 = ax.twinx()
    ax2.plot(NA_GRID, opt_density, "s--", color="tab:red", alpha=0.7)
    ax2.set_ylabel("optimal labeling density", color="tab:red")
    ax2.tick_params(axis="y", labelcolor="tab:red")
    ax.set_title("Best achievable usable-synapse yield rises with NA\n(red = density that achieves it)")
    ax.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig("figures/peak_yield_vs_na.png", dpi=120)

    # ---- Fig C: heatmap density x NA ----
    fig, ax = plt.subplots(figsize=(8.5, 6))
    pcm = ax.pcolormesh(DENSITIES, NA_GRID, total.T, cmap="magma", shading="auto")
    cs = ax.contour(DENSITIES, NA_GRID, total.T, levels=[1000, 2000, 3000, 4000, 5000],
                    colors="white", linewidths=0.8)
    ax.clabel(cs, fmt="%d", fontsize=8)
    ax.axhline(0.6, color="cyan", ls=":", lw=1)
    fig.colorbar(pcm, ax=ax, label="resolvable synapses (10 FOVs)")
    ax.set_xlabel("labeling density"); ax.set_ylabel("NA")
    ax.set_title("Usable-synapse yield over density x NA (white = iso-count contours)")
    fig.tight_layout(); fig.savefig("figures/yield_heatmap.png", dpi=120)

    print(f"peak yield {total.max():.0f} at density {opt_density[np.argmax(peak)]:.2f}, "
          f"NA {NA_GRID[np.argmax(peak)]:.2f}")
    print("peak vs NA:", {float(na): int(p) for na, p in zip(NA_GRID, peak)})
    print("-> figures/yield_vs_density_byNA.png, peak_yield_vs_na.png, yield_heatmap.png")


if __name__ == "__main__":
    main()
