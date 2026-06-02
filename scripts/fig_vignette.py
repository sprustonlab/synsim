"""10-axon vignette: ground truth vs imaged optical section at three NAs.

Picks a typical (median-density) neighborhood of axons, all labeled, and marks
each bouton resolvable (green) or confused (red) at the z0 optical section.

    python scripts/fig_vignette.py
"""
import numpy as np
from scipy.spatial import cKDTree
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from synsim import build_scene, footprint, psf_sigmas_um, PIX_UM
from synsim.metrics import AXIAL_DETECT, PURITY_THRESH

N_AXONS = 10
OUT = "figures/vignette.png"


def main():
    rng = np.random.default_rng(0)
    scene = build_scene(crop_um=20.0, rng=rng)
    data = [(sid, v, b) for sid, v, b in scene.axons if len(b)]

    cent = np.array([v[:, :2].mean(0) for _, v, _ in data])
    dists, idxs = cKDTree(cent).query(cent, k=N_AXONS)
    seed = int(np.argsort(dists[:, -1])[len(dists) // 2])     # typical neighborhood
    group = [data[i] for i in idxs[seed]]
    print(f"vignette axons: {[g[0] for g in group]}")

    allv = np.concatenate([v for _, v, _ in group])
    lo = allv[:, :2].min(0) - 0.7
    span = (allv[:, :2].max(0) + 0.7 - lo).max()
    bz = np.concatenate([b[:, 2] for _, _, b in group])
    z0 = float(bz[np.argmax([(np.abs(bz - z) <= 1.0).sum() for z in bz])])
    print(f"view {span:.1f} um, plane z0 = {z0:.1f} um")

    colors = plt.cm.tab10(np.linspace(0, 1, 10))
    fig, ax = plt.subplots(1, 4, figsize=(18, 5))

    for k, (sid, v, b) in enumerate(group):
        ax[0].scatter(v[:, 1] - lo[1], v[:, 0] - lo[0], s=0.5, color=colors[k])
        ax[0].scatter(b[:, 1] - lo[1], b[:, 0] - lo[0], s=55, facecolor=colors[k],
                      edgecolor="white", linewidth=1.2, zorder=5)
    ax[0].set_xlim(0, span); ax[0].set_ylim(0, span); ax[0].set_aspect("equal")
    ax[0].set_title("ground truth\n(dots=axon, ringed=boutons)", fontsize=10)
    ax[0].set_xlabel("um"); ax[0].set_ylabel("um")

    for col, na in zip(range(1, 4), (0.8, 0.6, 0.4)):
        _, s_z = psf_sigmas_um(na)
        fps = [footprint(v, lo, na, z0, span) for _, v, _ in group]
        total = np.sum(fps, axis=0)
        ax[col].imshow(total, origin="lower", cmap="gray",
                       extent=[0, span, 0, span])
        nres = ncon = 0
        for k, (sid, v, b) in enumerate(group):
            aw = np.exp(-((b[:, 2] - z0) ** 2) / (2 * s_z ** 2))
            for j in range(len(b)):
                if aw[j] <= AXIAL_DETECT:
                    continue
                gx = int((b[j, 0] - lo[0]) / PIX_UM)
                gy = int((b[j, 1] - lo[1]) / PIX_UM)
                if not (0 <= gx < total.shape[0] and 0 <= gy < total.shape[1]):
                    continue
                tv = total[gx, gy]; pv = fps[k][gx, gy]
                ok = tv > 0 and pv / tv >= PURITY_THRESH
                nres += ok; ncon += (not ok)
                ax[col].scatter(b[j, 1] - lo[1], b[j, 0] - lo[0], s=55, facecolor="none",
                                edgecolor=("lime" if ok else "red"), linewidth=1.4, zorder=5)
        tot = nres + ncon
        pct = f"{100*nres/tot:.0f}%" if tot else "n/a"
        ax[col].set_title(f"imaged NA {na}\n{nres}/{tot} in-plane resolvable ({pct})", fontsize=10)
        ax[col].set_xlabel("um")

    leg = [Line2D([0], [0], marker="o", color="w", markerfacecolor="none",
                  markeredgecolor="lime", markersize=10, label="resolvable (parent >70% of signal)"),
           Line2D([0], [0], marker="o", color="w", markerfacecolor="none",
                  markeredgecolor="red", markersize=10, label="confused (merged with neighbor)")]
    fig.legend(handles=leg, loc="lower center", ncol=2, fontsize=10, frameon=False)
    fig.suptitle("10 real LICONN axons: synapse assignability collapses as NA drops "
                 "(all 10 labeled)", fontsize=12)
    fig.tight_layout(rect=[0, 0.05, 1, 0.96])
    fig.savefig(OUT, dpi=110)
    print(f"-> {OUT}")


if __name__ == "__main__":
    main()
