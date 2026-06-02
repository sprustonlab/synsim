"""Two analysis figures over a density sweep, averaged across 10 FOVs.

Fig 1 (figures/na_for_50.png): for each labeling density (0.01-0.3, 15 points),
the NA needed to resolve 50% of imaged synapses, mean +/- std over 10 FOVs.

Fig 2 (figures/purity_dist.png): for each density, the distribution of per-bouton
purity at NA 0.6 (pooled over the 10 FOVs); each violin annotated with the number
of boutons with purity > 0.70 (= resolvable count).

    python scripts/fig_resolution_analysis.py
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
REGIONS = list(range(1, 11))          # 10 random sub-volumes
NA_FIG2 = 0.6
TARGET = 0.50


def crossing(na_grid, fracs, target):
    """Smallest NA at which resolvable fraction reaches `target` (lin. interp)."""
    g = np.array([a for a, f in zip(na_grid, fracs) if not np.isnan(f)])
    f = np.array([f for f in fracs if not np.isnan(f)])
    if len(f) < 2:
        return np.nan
    if f[0] >= target:
        return g[0]                   # already resolved at the lowest NA
    if f[-1] < target:
        return np.nan                 # not achievable within NA <= 1.0
    i = int(np.argmax(f >= target))
    x0, x1, y0, y1 = g[i - 1], g[i], f[i - 1], f[i]
    return x1 if y1 == y0 else x0 + (target - y0) * (x1 - x0) / (y1 - y0)


def _work(region):
    """One FOV: build scene, sweep density x NA, return NA-for-50% row + purity lists."""
    scene = build_scene(rng=np.random.default_rng(0),
                        params=replace(DEFAULTS, region=float(region)))
    u = np.random.default_rng(100 + region).random(len(scene.axons))
    na_row = np.full(len(DENSITIES), np.nan)
    purity = [[] for _ in DENSITIES]
    for di, dens in enumerate(DENSITIES):
        mask = u < dens
        fracs = []
        for na in NA_GRID:
            r = evaluate(scene, mask, replace(DEFAULTS, na=float(na)))
            fracs.append(r["resolvable_frac"] if r["scored"] > 0 else np.nan)
        na_row[di] = crossing(NA_GRID, fracs, TARGET)
        r06 = evaluate(scene, mask, replace(DEFAULTS, na=NA_FIG2))
        purity[di] = r06["purity"].tolist()
    return region, len(scene.axons), na_row, purity


def main():
    nproc = min(len(REGIONS), os.cpu_count() or 1)
    print(f"running {len(REGIONS)} FOVs on {nproc} cores...", flush=True)
    with Pool(nproc) as pool:
        results = pool.map(_work, REGIONS)

    na_needed = np.full((len(REGIONS), len(DENSITIES)), np.nan)
    purity_pool = [[] for _ in DENSITIES]
    for ri, (region, n_ax, na_row, purity) in enumerate(results):
        na_needed[ri] = na_row
        for di in range(len(DENSITIES)):
            purity_pool[di].extend(purity[di])
        print(f"region {region}: {n_ax} axons", flush=True)

    # ---- Fig 1: NA needed for 50% vs density ----
    mean = np.nanmean(na_needed, axis=0)
    std = np.nanstd(na_needed, axis=0)
    n_ok = np.sum(~np.isnan(na_needed), axis=0)
    fig, ax = plt.subplots(figsize=(8, 5))
    for ri in range(len(REGIONS)):
        ax.scatter(DENSITIES, na_needed[ri], s=10, color="0.8", zorder=1)
    ax.errorbar(DENSITIES, mean, yerr=std, fmt="o-", color="tab:blue", capsize=3, zorder=3)
    ax.set_xlabel("labeling density (fraction of axons)")
    ax.set_ylabel("NA needed to resolve 50% of imaged synapses")
    ax.set_ylim(0.3, 1.05)
    ax.set_title("NA for 50% resolvable vs density (mean +/- std, 10 FOVs)")
    ax.grid(alpha=0.3)
    for x, m, k in zip(DENSITIES, mean, n_ok):
        if k < len(REGIONS):
            ax.annotate(f"{k}/10", (x, 1.02), fontsize=7, ha="center", color="gray")
    fig.tight_layout(); fig.savefig("figures/na_for_50.png", dpi=120)

    # ---- Fig 2: purity distribution per density at NA 0.6 ----
    fig, ax = plt.subplots(figsize=(11, 5))
    data = [np.array(p) if len(p) else np.array([np.nan]) for p in purity_pool]
    pos = np.arange(len(DENSITIES))
    vp = ax.violinplot(data, positions=pos, showmedians=True, widths=0.8)
    for b in vp["bodies"]:
        b.set_facecolor("tab:green"); b.set_alpha(0.4)
    ax.axhline(0.70, color="red", ls="--", lw=1, label="purity 0.70 threshold")
    for i, p in enumerate(purity_pool):
        p = np.array(p)
        n_hi = int((p > 0.70).sum())
        ax.annotate(f"{n_hi}", (i, 1.04), fontsize=8, ha="center", color="darkgreen")
    ax.set_xticks(pos)
    ax.set_xticklabels([f"{d:.02f}" for d in DENSITIES], rotation=45, fontsize=8)
    ax.set_xlabel("labeling density")
    ax.set_ylabel("per-bouton purity")
    ax.set_ylim(0, 1.1)
    ax.set_title(f"Purity distribution at NA {NA_FIG2} vs density (10 FOVs pooled); "
                 f"numbers above = boutons with purity > 0.70")
    ax.legend(loc="lower left")
    fig.tight_layout(); fig.savefig("figures/purity_dist.png", dpi=120)

    # ---- printed summary ----
    print("\n density   NA_for_50 (mean+-std, n_ok)   #purity>0.70 (NA0.6, pooled)")
    for di, d in enumerate(DENSITIES):
        nhi = int((np.array(purity_pool[di]) > 0.70).sum())
        ntot = len(purity_pool[di])
        print(f" {d:5.3f}    {mean[di]:.3f} +- {std[di]:.3f}  ({n_ok[di]}/10)      "
              f"{nhi}/{ntot}")
    print("-> figures/na_for_50.png, figures/purity_dist.png")


if __name__ == "__main__":
    main()
