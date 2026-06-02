"""Redo the 10-axon confusion figure with the TRUE 3D imaging core, and test
whether the old 2D `footprint` matches a real 3D-convolution optical section.

Outputs:
  spike_equiv.png       - old 2D footprint vs 3D z0-section (NA 0.6) + difference
  spike_demo10_3d.png   - 10-axon vignette imaged with the 3D engine (top section
                          + side projection), boutons resolvable/confused.
Run: .venv-liconn/bin/python spike_demo10_3d.py
"""
import numpy as np
from scipy.spatial import cKDTree
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from spike_sim import (load_crop_axons, skeleton_and_boutons, psf_sigmas_um,
                       footprint, PIX_UM, AXIAL_DETECT, PURITY_THRESH)
from imaging import voxelize, blur3d, section_index

GRID = 0.10
N_AXONS = 10


def select_group():
    axons, lo, hi = load_crop_axons()
    data = []
    for sid, v in axons:
        L, b = skeleton_and_boutons(v)
        if len(b):
            data.append((sid, v, b))
    cent = np.array([v[:, :2].mean(0) for _, v, _ in data])
    tree = cKDTree(cent)
    dists, idxs = tree.query(cent, k=N_AXONS)
    seed = int(np.argsort(dists[:, -1])[len(dists) // 2])  # median neighborhood
    group = [data[i] for i in idxs[seed]]
    return group


def build_grid(group, pad=2.0):
    allv = np.concatenate([v for _, v, _ in group])
    origin = allv.min(0) - pad
    shape = tuple((np.ceil((allv.max(0) + pad - origin) / GRID)).astype(int))
    return origin, shape


def per_axon_3d(group, origin, shape, na):
    """Return list of blurred 3D volumes (one per axon) + their sum."""
    vols = []
    total = np.zeros(shape, np.float32)
    for sid, v, b in group:
        bl = blur3d(voxelize(v, origin, shape, GRID), na, GRID)
        vols.append(bl)
        total += bl
    return vols, total


def equivalence_test(group, origin, shape):
    """Compare old 2D footprint (point cloud, per-vertex axial weight) vs the
    3D-convolution optical section at the same plane, NA 0.6."""
    na = 0.6
    bz = np.concatenate([b[:, 2] for _, _, b in group])
    z0 = float(bz[np.argmax([(np.abs(bz - z) <= 1.0).sum() for z in bz])])
    k0 = section_index(z0, origin, GRID)

    # 3D engine: sum of per-axon blurred volumes, take z0 slice
    _, total3d = per_axon_3d(group, origin, shape, na)
    sec3d = total3d[:, :, k0]                       # (x, y) on GRID

    # old 2D footprint engine, on the SAME lateral window/pixel as `footprint`
    lo2d = origin[:2]
    npx2 = int(np.ceil((shape[0] * GRID) / PIX_UM))
    # footprint() uses module CROP grid; recompute a matching 2D image here
    s_xy, s_z = psf_sigmas_um(na)
    from scipy import ndimage
    img2d = np.zeros((shape[0], shape[1]), np.float32)
    for sid, v, b in group:
        w = np.exp(-((v[:, 2] - z0) ** 2) / (2 * s_z ** 2))
        gx = ((v[:, 0] - origin[0]) / GRID).astype(int)
        gy = ((v[:, 1] - origin[1]) / GRID).astype(int)
        ok = (gx >= 0) & (gx < shape[0]) & (gy >= 0) & (gy < shape[1])
        np.add.at(img2d, (gx[ok], gy[ok]), w[ok])
    img2d = ndimage.gaussian_filter(img2d, s_xy / GRID)

    # normalize each to its own max for shape comparison
    a = sec3d / (sec3d.max() + 1e-9)
    bimg = img2d / (img2d.max() + 1e-9)
    mask = (a > 0.02) | (bimg > 0.02)
    corr = np.corrcoef(a[mask], bimg[mask])[0, 1]
    rel = np.abs(a - bimg)[mask].mean()
    print(f"[equiv NA0.6] Pearson r = {corr:.3f} | mean|diff| (norm) = {rel:.3f}")

    fig, ax = plt.subplots(1, 3, figsize=(15, 5))
    ext = [0, shape[1] * GRID, 0, shape[0] * GRID]
    ax[0].imshow(bimg.T, origin="lower", cmap="magma", extent=ext)
    ax[0].set_title("OLD 2D footprint\n(surface points, per-vertex axial weight)")
    ax[1].imshow(a.T, origin="lower", cmap="magma", extent=ext)
    ax[1].set_title("3D conv -> z0 optical section\n(filled tube, true 3D PSF)")
    im = ax[2].imshow((a - bimg).T, origin="lower", cmap="RdBu", vmin=-1, vmax=1, extent=ext)
    ax[2].set_title(f"difference (3D - 2D)\nPearson r={corr:.3f}, mean|diff|={rel:.3f}")
    for a_ in ax:
        a_.set_xlabel("y (um)"); a_.set_ylabel("x (um)")
    fig.colorbar(im, ax=ax[2], fraction=0.046)
    fig.suptitle("Do the two imaging models agree? (same 10 axons, same plane, NA 0.6)",
                 fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig("spike_equiv.png", dpi=110)
    print("-> spike_equiv.png")
    return z0, k0


def redo_demo10(group, origin, shape, z0, k0):
    cols = ["ground truth"] + [f"NA {na}" for na in (0.8, 0.6, 0.4)]
    fig, ax = plt.subplots(2, len(cols), figsize=(4.4 * len(cols), 8),
                           gridspec_kw={"height_ratios": [shape[1], shape[2]]})
    nx, ny, nz = shape
    ext_top = [0, nx * GRID, 0, ny * GRID]
    ext_side = [0, nx * GRID, 0, nz * GRID]

    # ground truth (no PSF): MIP of voxelized tubes
    gt = np.zeros(shape, np.float32)
    for sid, v, b in group:
        gt = np.maximum(gt, voxelize(v, origin, shape, GRID))
    ax[0, 0].imshow(gt.max(2).T, origin="lower", cmap="magma", aspect="equal", extent=ext_top)
    ax[1, 0].imshow(gt.max(1).T, origin="lower", cmap="magma", aspect="equal", extent=ext_side)
    ax[0, 0].set_title("ground truth (MIP)\nringed = boutons", fontsize=10)
    for sid, v, b in group:
        ax[0, 0].scatter(b[:, 0] - origin[0], b[:, 1] - origin[1], s=40,
                         facecolor="none", edgecolor="cyan", lw=1.0, zorder=5)

    for c, na in enumerate((0.8, 0.6, 0.4), start=1):
        s_xy, s_z = psf_sigmas_um(na)
        vols, total = per_axon_3d(group, origin, shape, na)
        sec = total[:, :, k0]                          # z0 optical section (top)
        ax[0, c].imshow(sec.T, origin="lower", cmap="gray", aspect="equal", extent=ext_top)
        ax[1, c].imshow(total.max(1).T, origin="lower", cmap="gray", aspect="equal", extent=ext_side)
        ax[1, c].axhline(z0 - origin[2], color="cyan", lw=0.8, ls=":")
        nres = ncon = 0
        for k, (sid, v, b) in enumerate(group):
            aw = np.exp(-((b[:, 2] - z0) ** 2) / (2 * s_z ** 2))
            psec = vols[k][:, :, k0]
            for j in range(len(b)):
                if aw[j] <= AXIAL_DETECT:
                    continue
                gx = int((b[j, 0] - origin[0]) / GRID)
                gy = int((b[j, 1] - origin[1]) / GRID)
                if not (0 <= gx < nx and 0 <= gy < ny):
                    continue
                tv = sec[gx, gy]; pv = psec[gx, gy]
                ok = tv > 0 and pv / tv >= PURITY_THRESH
                nres += ok; ncon += (not ok)
                ax[0, c].scatter(b[j, 0] - origin[0], b[j, 1] - origin[1], s=40,
                                 facecolor="none",
                                 edgecolor=("lime" if ok else "red"), lw=1.3, zorder=5)
        tot = nres + ncon
        pct = f"{100*nres/tot:.0f}%" if tot else "n/a"
        ax[0, c].set_title(f"NA {na}: z0 section\n{nres}/{tot} resolvable ({pct})", fontsize=10)
    for c in range(len(cols)):
        ax[0, c].set_xlabel("x (um)"); ax[0, c].set_ylabel("y (um)")
        ax[1, c].set_xlabel("x (um)"); ax[1, c].set_ylabel("z (um)")

    leg = [Line2D([0], [0], marker="o", color="w", markerfacecolor="none",
                  markeredgecolor="lime", markersize=10, label="resolvable (parent >70% at z0)"),
           Line2D([0], [0], marker="o", color="w", markerfacecolor="none",
                  markeredgecolor="red", markersize=10, label="confused (merged at z0)")]
    fig.legend(handles=leg, loc="lower center", ncol=2, fontsize=10, frameon=False)
    fig.suptitle("10 real LICONN axons imaged with the 3D engine: top = z0 optical section, "
                 "bottom = side projection (z-smear). All 10 labeled.", fontsize=12)
    fig.tight_layout(rect=[0, 0.04, 1, 0.95])
    fig.savefig("spike_demo10_3d.png", dpi=110)
    print("-> spike_demo10_3d.png")


def main():
    group = select_group()
    print(f"ids {[g[0] for g in group]}")
    origin, shape = build_grid(group)
    print(f"shared grid {shape} @ {GRID} um")
    z0, k0 = equivalence_test(group, origin, shape)
    print(f"plane z0={z0:.1f} um (k0={k0})")
    redo_demo10(group, origin, shape, z0, k0)


if __name__ == "__main__":
    main()
