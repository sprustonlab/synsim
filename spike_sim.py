"""SPIKE: sparse-axon imaging confusion, on real LICONN axon geometry.

Proof-of-concept only. Validates the full core before any build:
  load axons in a crop -> 3D skeletonize -> length -> place boutons at ~4.3 um
  spacing -> image one plane at a given NA -> score each bouton's resolvability
  (parent-axon signal fraction at its pixel) across labeling density x NA.

Run: .venv-liconn/bin/python spike_sim.py
"""
import glob
import time

import numpy as np
from scipy import ndimage
from skimage.morphology import skeletonize
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RES_NM = np.array([18.0, 18.0, 24.0])            # native voxel size
UM = RES_NM / 1000.0                              # voxel -> um per axis
CROP_UM = 20.0                                    # lateral crop side
PIX_UM = 0.15                                     # sim image pixel
SKEL_UM = 0.25                                    # voxel size for skeletonization
BOUTON_SPACING_UM = 4.3                           # calibrated: 0.95/4.084 = 0.233/um
PURITY_THRESH = 0.70
AXIAL_DETECT = 0.5                                 # bouton "in-plane" if axial weight > this
RNG = np.random.default_rng(0)


def psf_sigmas_um(na):
    fwhm_xy = 0.66 * (0.6 / na)                    # lateral ~ 1/NA
    fwhm_z = 4.09 * (0.6 / na) ** 2                # axial  ~ 1/NA^2
    return fwhm_xy / 2.3548, fwhm_z / 2.3548


def load_crop_axons():
    """Return list of (id, verts_um Nx3) for axons with vertices in a central crop."""
    files = sorted(glob.glob("data/axons/chunk_*.npz"))
    # first pass: global lateral centroid to center the crop
    sx = sy = n = 0
    for f in files:
        v = np.load(f)["verts"]
        if v.shape[0]:
            sx += v[:, 0].sum() * UM[0]
            sy += v[:, 1].sum() * UM[1]
            n += v.shape[0]
    cx, cy = sx / n, sy / n
    lo = np.array([cx - CROP_UM / 2, cy - CROP_UM / 2])
    hi = lo + CROP_UM
    axons = []
    for f in files:
        d = np.load(f)
        verts, offs, ids = d["verts"], d["offsets"], d["ids"]
        for k, sid in enumerate(ids):
            v = verts[offs[k]:offs[k + 1]].astype(np.float64) * UM   # -> um
            inx = (v[:, 0] >= lo[0]) & (v[:, 0] < hi[0]) & (v[:, 1] >= lo[1]) & (v[:, 1] < hi[1])
            if inx.sum() >= 50:                    # enough of the axon is in-crop
                axons.append((int(sid), v[inx]))
    return axons, lo, hi


def skeleton_and_boutons(v_um):
    """3D skeletonize an axon's points; return (length_um, bouton_xyz Nx3)."""
    mn = v_um.min(0)
    idx = np.floor((v_um - mn) / SKEL_UM).astype(int)
    shp = idx.max(0) + 3
    vol = np.zeros(shp, bool)
    vol[idx[:, 0] + 1, idx[:, 1] + 1, idx[:, 2] + 1] = True
    vol = ndimage.binary_dilation(vol, iterations=1)
    vol = ndimage.binary_fill_holes(vol)
    skel = skeletonize(vol)
    sp = np.argwhere(skel)
    if sp.shape[0] == 0:
        return 0.0, np.zeros((0, 3))
    length = sp.shape[0] * SKEL_UM
    skel_um = (sp - 1) * SKEL_UM + mn
    nb = max(1, int(round(length / BOUTON_SPACING_UM)))
    pick = RNG.choice(sp.shape[0], size=min(nb, sp.shape[0]), replace=False)
    return length, skel_um[pick]


def footprint(v_um, lo, na, z0):
    """Axial-weighted, laterally-blurred 2D image of one axon in the crop grid."""
    s_xy, s_z = psf_sigmas_um(na)
    w = np.exp(-((v_um[:, 2] - z0) ** 2) / (2 * s_z ** 2))
    gx = ((v_um[:, 0] - lo[0]) / PIX_UM).astype(int)
    gy = ((v_um[:, 1] - lo[1]) / PIX_UM).astype(int)
    npx = int(CROP_UM / PIX_UM)
    ok = (gx >= 0) & (gx < npx) & (gy >= 0) & (gy < npx)
    img = np.zeros((npx, npx))
    np.add.at(img, (gx[ok], gy[ok]), w[ok])
    return ndimage.gaussian_filter(img, s_xy / PIX_UM)


def main():
    t0 = time.time()
    axons, lo, hi = load_crop_axons()
    print(f"crop {CROP_UM}x{CROP_UM} um: {len(axons)} axons  ({time.time()-t0:.1f}s)")

    # skeletons + boutons
    tot_len = 0.0
    data = []  # (id, verts_um, boutons_xyz)
    for sid, v in axons:
        L, b = skeleton_and_boutons(v)
        tot_len += L
        data.append((sid, v, b))
    n_b = sum(len(b) for _, _, b in data)
    print(f"total axon length {tot_len:.0f} um | boutons {n_b} | "
          f"synapses/um = {n_b/tot_len:.3f} (target ~0.233)  ({time.time()-t0:.1f}s)")

    z0 = np.median(np.concatenate([b[:, 2] for _, _, b in data if len(b)]))
    print(f"imaging plane z0 = {z0:.1f} um\n")

    print(f"{'NA':>4} {'dens':>5} {'lab.ax':>7} {'det.bout':>9} "
          f"{'resolv':>7} {'resolv%':>8}")
    results = {}
    for na in (0.4, 0.6, 0.8):
        s_xy, s_z = psf_sigmas_um(na)
        for dens in (0.05, 0.20, 0.60):
            lab = RNG.random(len(data)) < dens
            total = np.zeros((int(CROP_UM / PIX_UM),) * 2)
            fps = {}
            for i, (sid, v, b) in enumerate(data):
                if lab[i]:
                    fp = footprint(v, lo, na, z0)
                    fps[i] = fp
                    total += fp
            det = resolv = 0
            for i, (sid, v, b) in enumerate(data):
                if not lab[i] or len(b) == 0:
                    continue
                aw = np.exp(-((b[:, 2] - z0) ** 2) / (2 * s_z ** 2))
                gx = ((b[:, 0] - lo[0]) / PIX_UM).astype(int)
                gy = ((b[:, 1] - lo[1]) / PIX_UM).astype(int)
                npx = total.shape[0]
                for j in range(len(b)):
                    if aw[j] <= AXIAL_DETECT or not (0 <= gx[j] < npx and 0 <= gy[j] < npx):
                        continue
                    det += 1
                    tv = total[gx[j], gy[j]]
                    pv = fps[i][gx[j], gy[j]]
                    if tv > 0 and pv / tv >= PURITY_THRESH:
                        resolv += 1
            frac = resolv / det if det else float("nan")
            results[(na, dens)] = (lab.sum(), det, resolv, frac)
            print(f"{na:>4.1f} {dens:>5.2f} {lab.sum():>7d} {det:>9d} "
                  f"{resolv:>7d} {100*frac:>7.1f}%")

    # demo PNG: two nearby axons, ground truth vs imaged at NA 0.6
    make_demo_png(data, lo, z0)
    print(f"\nDONE ({time.time()-t0:.1f}s) -> spike_demo.png")


def make_demo_png(data, lo, z0):
    cand = [(i, v, b) for i, (sid, v, b) in enumerate(data) if len(b)]
    if len(cand) < 2:
        return
    cents = {i: v[:, :2].mean(0) for i, v, b in cand}
    best = None
    items = list(cents.items())
    for a in range(len(items)):
        for c in range(a + 1, len(items)):
            d = np.hypot(*(items[a][1] - items[c][1]))
            if best is None or d < best[0]:
                best = (d, items[a][0], items[c][0])
    _, i1, i2 = best
    v1 = data[i1][1]; v2 = data[i2][1]
    b1 = data[i1][2]; b2 = data[i2][2]
    fig, ax = plt.subplots(1, 4, figsize=(17, 4.3))
    for k, na in enumerate([None, 0.8, 0.6, 0.4]):
        if na is None:
            ax[k].scatter((v1[:, 1]-lo[1]), (v1[:, 0]-lo[0]), s=0.4, c="tab:blue")
            ax[k].scatter((v2[:, 1]-lo[1]), (v2[:, 0]-lo[0]), s=0.4, c="tab:red")
            ax[k].set_title("ground truth (2 axons)")
        else:
            f1 = footprint(v1, lo, na, z0); f2 = footprint(v2, lo, na, z0)
            rgb = np.zeros(f1.shape + (3,))
            rgb[..., 2] = f1 / (f1.max() + 1e-9)
            rgb[..., 0] = f2 / (f2.max() + 1e-9)
            ax[k].imshow(np.clip(rgb, 0, 1), origin="lower")
            for b, col in [(b1, "cyan"), (b2, "yellow")]:
                gx = (b[:, 0]-lo[0])/PIX_UM; gy = (b[:, 1]-lo[1])/PIX_UM
                ax[k].scatter(gy, gx, s=14, edgecolor=col, facecolor="none", lw=0.8)
            ax[k].set_title(f"imaged NA {na}")
        ax[k].set_xticks([]); ax[k].set_yticks([])
    fig.suptitle("Real LICONN axons: as NA drops, axial+lateral PSF merges neighbors")
    fig.tight_layout()
    fig.savefig("spike_demo.png", dpi=110)


if __name__ == "__main__":
    main()
