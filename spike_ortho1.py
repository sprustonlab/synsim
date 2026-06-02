"""Single-axon orthogonal-view sanity check: top (xy) + side (xz) at several NAs.

Validates that the anisotropic 3D PSF is applied correctly: lateral blur stays
tight while axial blur is ~4 um at NA 0.6 (sigma_z ~ 6x sigma_xy). True 3D
voxelize -> anisotropic 3D Gaussian -> max-intensity projections.

Run: .venv-liconn/bin/python spike_ortho1.py
"""
import glob

import numpy as np
from scipy import ndimage
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from spike_sim import psf_sigmas_um

RES_NM = np.array([18.0, 18.0, 24.0])
UM = RES_NM / 1000.0
GRID = 0.10  # isotropic sim voxel (um) -> PSF anisotropy comes only from the PSF
NAS = (0.8, 0.6, 0.4)


def pick_axon():
    """Choose one axon with large z-extent (most informative side view)."""
    best = None
    for f in sorted(glob.glob("data/axons/chunk_*.npz")):
        d = np.load(f)
        verts, offs, ids = d["verts"], d["offsets"], d["ids"]
        for k, sid in enumerate(ids):
            v = verts[offs[k]:offs[k + 1]].astype(np.float64) * UM
            if v.shape[0] < 800:
                continue
            ext = v.max(0) - v.min(0)
            score = ext[2] * v.shape[0]            # z-extent * size
            if best is None or score > best[0]:
                best = (score, int(sid), v)
    return best[1], best[2]


def voxelize(v_um):
    mn = v_um.min(0)
    idx = np.floor((v_um - mn) / GRID).astype(int)
    shp = idx.max(0) + 3
    vol = np.zeros(shp, bool)
    vol[idx[:, 0] + 1, idx[:, 1] + 1, idx[:, 2] + 1] = True
    vol = ndimage.binary_dilation(vol, iterations=1)
    vol = ndimage.binary_fill_holes(vol)
    return vol.astype(np.float32), mn


def main():
    sid, v = pick_axon()
    ext = v.max(0) - v.min(0)
    print(f"axon {sid}: {v.shape[0]} verts, extent (um) "
          f"x={ext[0]:.1f} y={ext[1]:.1f} z={ext[2]:.1f}")
    vol, mn = voxelize(v)
    nx, ny, nz = vol.shape
    print(f"grid {nx}x{ny}x{nz} @ {GRID} um")

    cols = ["ground truth"] + [f"NA {na}" for na in NAS]
    fig, ax = plt.subplots(2, len(cols), figsize=(4.2 * len(cols), 7),
                           gridspec_kw={"height_ratios": [ny, nz]})

    def show(col, vol3d, title):
        top = vol3d.max(axis=2).T          # (y, x)
        side = vol3d.max(axis=1).T         # (z, x)
        ax[0, col].imshow(top, origin="lower", cmap="magma", aspect="equal",
                          extent=[0, nx * GRID, 0, ny * GRID])
        ax[1, col].imshow(side, origin="lower", cmap="magma", aspect="equal",
                          extent=[0, nx * GRID, 0, nz * GRID])
        ax[0, col].set_title(title, fontsize=11)
        ax[0, col].set_ylabel("y (um)"); ax[1, col].set_ylabel("z (um)")
        ax[1, col].set_xlabel("x (um)")

    show(0, vol, "ground truth")
    for c, na in enumerate(NAS, start=1):
        s_xy, s_z = psf_sigmas_um(na)
        blur = ndimage.gaussian_filter(vol, (s_xy / GRID, s_xy / GRID, s_z / GRID))
        show(c, blur, f"NA {na}\nsigma_xy={s_xy:.2f}  sigma_z={s_z:.2f} um")

    fig.text(0.5, 0.965, f"Axon {sid}: top (xy) vs side (xz). Lateral blur stays tight; "
             f"axial blur grows ~1/NA^2 (sigma_z ~ {psf_sigmas_um(0.6)[1]:.1f} um at NA 0.6)",
             ha="center", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig("spike_ortho1.png", dpi=110)
    print("-> spike_ortho1.png")


if __name__ == "__main__":
    main()
