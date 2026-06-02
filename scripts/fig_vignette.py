"""10-axon vignette: ground truth vs imaged section at three NAs.

Picks a typical neighborhood of axons (all labeled), scores every bouton with the
gridless 3D metric, and marks each resolvable (green) or confused (red). The
grayscale background is a display optical section (footprint at the FOV mid-plane);
the green/red verdicts come from the real 3D scoring, not the 2D display.

    python scripts/fig_vignette.py
"""
import numpy as np
from scipy.spatial import cKDTree
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from synsim import build_scene, footprint, evaluate, SimParams, PIX_UM
from synsim.scene import Scene

N_AXONS = 10
OUT = "figures/vignette.png"


def main():
    rng = np.random.default_rng(0)
    scene = build_scene(rng=rng)
    data = [(sid, v, b) for sid, v, b in scene.axons if len(b)]

    cent = np.array([b[:, :2].mean(0) for _, _, b in data])
    dists, idxs = cKDTree(cent).query(cent, k=N_AXONS)
    seed = int(np.argsort(dists[:, -1])[len(dists) // 2])     # typical neighborhood
    group = [data[i] for i in idxs[seed]]
    print(f"vignette axons: {[g[0] for g in group]}")

    allv = np.concatenate([v for _, v, _ in group])
    lo = allv.min(0) - 0.7
    fov = (allv.max(0) + 0.7) - lo
    span = float(fov[:2].max())
    z0 = float(lo[2] + fov[2] / 2)
    sub = Scene(group, lo, fov, total_length_um=1.0)
    mask = np.ones(len(group), bool)

    colors = plt.cm.tab10(np.linspace(0, 1, 10))
    fig, ax = plt.subplots(1, 4, figsize=(18, 5))

    for k, (sid, v, b) in enumerate(group):
        ax[0].scatter(v[:, 1] - lo[1], v[:, 0] - lo[0], s=0.6, color=colors[k])
        ax[0].scatter(b[:, 1] - lo[1], b[:, 0] - lo[0], s=55, facecolor=colors[k],
                      edgecolor="white", linewidth=1.2, zorder=5)
    ax[0].set_xlim(0, span); ax[0].set_ylim(0, span); ax[0].set_aspect("equal")
    ax[0].set_title("ground truth\n(dots=axon, ringed=boutons)", fontsize=10)
    ax[0].set_xlabel("um"); ax[0].set_ylabel("um")

    for col, na in zip(range(1, 4), (0.8, 0.6, 0.4)):
        params = SimParams(na=na, edge_xy_um=0.0, edge_z_um=0.0)
        r = evaluate(sub, mask, params)
        total = np.sum([footprint(v, lo[:2], na, z0, span) for _, v, _ in group], axis=0)
        ax[col].imshow(total, origin="lower", cmap="gray", extent=[0, span, 0, span])
        pos, res = r["positions"], r["resolved"]
        for good, c in [(res, "lime"), (~res, "red")]:
            ax[col].scatter(pos[good, 1] - lo[1], pos[good, 0] - lo[0], s=55,
                            facecolor="none", edgecolor=c, linewidth=1.4, zorder=5)
        frac = 100 * r["resolvable_frac"] if r["scored"] else float("nan")
        ax[col].set_title(f"imaged NA {na}\n{r['resolvable']}/{r['scored']} resolvable "
                          f"({frac:.0f}%)", fontsize=10)
        ax[col].set_xlabel("um")

    leg = [Line2D([0], [0], marker="o", color="w", markerfacecolor="none",
                  markeredgecolor="lime", markersize=10, label="resolvable (parent >70% of 3D signal)"),
           Line2D([0], [0], marker="o", color="w", markerfacecolor="none",
                  markeredgecolor="red", markersize=10, label="confused (merged with neighbor)")]
    fig.legend(handles=leg, loc="lower center", ncol=2, fontsize=10, frameon=False)
    fig.suptitle("10 real LICONN axons: synapse assignability collapses as NA drops "
                 "(all 10 labeled; verdicts from gridless 3D scoring)", fontsize=12)
    fig.tight_layout(rect=[0, 0.05, 1, 0.96])
    fig.savefig(OUT, dpi=110)
    print(f"-> {OUT}")


if __name__ == "__main__":
    main()
