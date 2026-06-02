"""Labeled 10-axon vignette for the spike: ground truth vs imaged at 3 NAs.

Reuses the spike's geometry + PSF helpers. Shows, in a tight view of 10
neighboring real LICONN axons:
  - ground truth: each axon a distinct color; boutons = white-ringed dots
  - imaged (NA 0.8/0.6/0.4): grayscale imaged plane; each bouton ringed
    GREEN if resolvable (parent axon > 70% of signal at its pixel) or
    RED if confused (merged with a neighbor).
Run: .venv-liconn/bin/python spike_demo10.py
"""
import numpy as np
from scipy import ndimage
from scipy.spatial import cKDTree
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from spike_sim import (load_crop_axons, skeleton_and_boutons, psf_sigmas_um,
                       PIX_UM, AXIAL_DETECT, PURITY_THRESH)

N_AXONS = 10


def local_footprint(v_um, origin, npx, na, z0):
    s_xy, s_z = psf_sigmas_um(na)
    w = np.exp(-((v_um[:, 2] - z0) ** 2) / (2 * s_z ** 2))
    gx = ((v_um[:, 0] - origin[0]) / PIX_UM).astype(int)
    gy = ((v_um[:, 1] - origin[1]) / PIX_UM).astype(int)
    ok = (gx >= 0) & (gx < npx) & (gy >= 0) & (gy < npx)
    img = np.zeros((npx, npx))
    np.add.at(img, (gx[ok], gy[ok]), w[ok])
    return ndimage.gaussian_filter(img, s_xy / PIX_UM)


def main():
    axons, lo, hi = load_crop_axons()
    data = []
    for sid, v in axons:
        L, b = skeleton_and_boutons(v)
        if len(b):
            data.append((sid, v, b))

    # pick the 10 tightest-clustered axons (so overlap/merging is visible)
    cent = np.array([v[:, :2].mean(0) for _, v, _ in data])
    tree = cKDTree(cent)
    dists, idxs = tree.query(cent, k=N_AXONS)
    # typical (median-density) neighborhood, not the tightest sub-PSF bundle
    seed = int(np.argsort(dists[:, -1])[len(dists) // 2])
    sel = idxs[seed]
    group = [data[i] for i in sel]
    print(f"selected {len(group)} axons; ids {[g[0] for g in group]}")

    allv = np.concatenate([v for _, v, _ in group])
    mn = allv[:, :2].min(0) - 0.7
    mx = allv[:, :2].max(0) + 0.7
    npx = int(np.ceil((mx - mn).max() / PIX_UM))
    origin = mn
    # plane = bouton-z with the most neighbors within +/-1 um (a populated plane)
    bz = np.concatenate([b[:, 2] for _, _, b in group])
    z0 = float(bz[np.argmax([(np.abs(bz - z) <= 1.0).sum() for z in bz])])
    print(f"view {npx*PIX_UM:.1f} um, plane z0={z0:.1f} um")

    colors = plt.cm.tab10(np.linspace(0, 1, 10))
    fig, ax = plt.subplots(1, 4, figsize=(18, 5))

    # ground truth
    for k, (sid, v, b) in enumerate(group):
        ax[0].scatter(v[:, 1] - origin[1], v[:, 0] - origin[0], s=0.5, color=colors[k])
        ax[0].scatter(b[:, 1] - origin[1], b[:, 0] - origin[0], s=55,
                      facecolor=colors[k], edgecolor="white", linewidth=1.2, zorder=5)
    ax[0].set_xlim(0, (mx - mn)[1]); ax[0].set_ylim(0, (mx - mn)[0])
    ax[0].set_aspect("equal")
    ax[0].set_title("ground truth\n(dots=axon membrane, ringed=boutons/synapses)", fontsize=10)
    ax[0].set_xlabel("um"); ax[0].set_ylabel("um")

    for col, na in zip(range(1, 4), (0.8, 0.6, 0.4)):
        s_xy, s_z = psf_sigmas_um(na)
        fps = [local_footprint(v, origin, npx, na, z0) for _, v, _ in group]
        total = np.sum(fps, axis=0)
        ax[col].imshow(total, origin="lower", cmap="gray",
                       extent=[0, (mx-mn)[1], 0, (mx-mn)[0]])
        nres = ncon = 0
        for k, (sid, v, b) in enumerate(group):
            aw = np.exp(-((b[:, 2] - z0) ** 2) / (2 * s_z ** 2))
            for j in range(len(b)):
                if aw[j] <= AXIAL_DETECT:
                    continue
                gx = int((b[j, 0] - origin[0]) / PIX_UM)
                gy = int((b[j, 1] - origin[1]) / PIX_UM)
                if not (0 <= gx < npx and 0 <= gy < npx):
                    continue
                tv = total[gx, gy]; pv = fps[k][gx, gy]
                ok = tv > 0 and pv / tv >= PURITY_THRESH
                nres += ok; ncon += (not ok)
                ax[col].scatter(b[j, 1] - origin[1], b[j, 0] - origin[0], s=55,
                                facecolor="none",
                                edgecolor=("lime" if ok else "red"), linewidth=1.4, zorder=5)
        tot = nres + ncon
        pct = f"{100*nres/tot:.0f}%" if tot else "n/a"
        ax[col].set_title(f"imaged NA {na}\n{nres}/{tot} in-plane boutons resolvable "
                          f"({pct})", fontsize=10)
        ax[col].set_xlabel("um")

    leg = [Line2D([0], [0], marker="o", color="w", markerfacecolor="none",
                  markeredgecolor="lime", markersize=10, label="resolvable bouton (parent >70% of signal)"),
           Line2D([0], [0], marker="o", color="w", markerfacecolor="none",
                  markeredgecolor="red", markersize=10, label="confused bouton (merged with neighbor)")]
    fig.legend(handles=leg, loc="lower center", ncol=2, fontsize=10, frameon=False)
    fig.suptitle("10 neighboring real LICONN axons: synapse assignability collapses as NA drops "
                 "(all 10 labeled in this vignette)", fontsize=12)
    fig.tight_layout(rect=[0, 0.05, 1, 0.96])
    fig.savefig("spike_demo10.png", dpi=110)
    print("-> spike_demo10.png")


if __name__ == "__main__":
    main()
