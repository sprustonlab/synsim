"""Show bouton placement on a single axon: 9 independent draws in a 3x3 grid.

Picks one full axon, then calls the bouton placer with 9 different RNG seeds.
The skeleton (centerline) and bouton count are fixed; only which skeleton points
are chosen varies, so the grid shows the stochasticity of placement.

    python scripts/fig_bouton_placements.py
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from synsim.data import iter_objects
from synsim.boutons import skeleton_and_boutons, skeleton_points, BOUTON_SPACING_UM

OUT = "figures/bouton_placements.png"


def pick_axon():
    """A full axon of moderate length (35-60 um) and modest lateral spread."""
    for sid, v in iter_objects("axons"):
        if not (2000 <= v.shape[0] <= 9000):
            continue
        ext = v.max(0) - v.min(0)
        if ext[:2].max() > 30:
            continue
        L, _ = skeleton_and_boutons(v, np.random.default_rng(0))
        if 35 <= L <= 60:
            return sid, v, L
    raise RuntimeError("no suitable axon found")


def main():
    sid, v, L = pick_axon()
    skel = skeleton_points(v)
    mn = v[:, :2].min(0)
    nb = max(1, int(round(L / BOUTON_SPACING_UM)))
    print(f"axon {sid}: length {L:.1f} um -> {nb} boutons/placement "
          f"(spacing {BOUTON_SPACING_UM} um); {len(skel)} skeleton points")

    fig, axes = plt.subplots(3, 3, figsize=(12, 12), sharex=True, sharey=True)
    for i, ax in enumerate(axes.flat):
        _, b = skeleton_and_boutons(v, np.random.default_rng(i))
        ax.scatter(v[:, 0] - mn[0], v[:, 1] - mn[1], s=0.6, color="0.8", zorder=1)
        ax.scatter(skel[:, 0] - mn[0], skel[:, 1] - mn[1], s=1.2, color="tab:blue", zorder=2)
        ax.scatter(b[:, 0] - mn[0], b[:, 1] - mn[1], s=90, facecolor="none",
                   edgecolor="red", linewidth=1.6, zorder=3)
        ax.set_title(f"placement seed {i}: {len(b)} boutons", fontsize=10)
        ax.set_aspect("equal")
    for ax in axes[-1, :]:
        ax.set_xlabel("x (um)")
    for ax in axes[:, 0]:
        ax.set_ylabel("y (um)")
    fig.suptitle(f"Axon {sid} (top view): boutons (red) sampled along skeleton (blue), "
                 f"gray=axon membrane. 9 independent placements.", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(OUT, dpi=110)
    print(f"-> {OUT}")


if __name__ == "__main__":
    main()
