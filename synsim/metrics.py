"""Resolvability metrics: image a labeled subset, score each bouton.

A bouton is 'detected' (in the optical section) if its parent axon's axial PSF
weight at the plane exceeds AXIAL_DETECT, and 'resolvable' if the parent axon
contributes more than PURITY_THRESH of the signal at the bouton's pixel.
"""
import numpy as np

from .imaging import footprint, PIX_UM
from .psf import psf_sigmas_um

AXIAL_DETECT = 0.5
PURITY_THRESH = 0.70


def image_labeled(scene, labeled_mask, na, z0, pix_um=PIX_UM):
    """Return (total_image, {axon_index: footprint}) over labeled axons."""
    npx = int(round(scene.crop_um / pix_um))
    total = np.zeros((npx, npx))
    fps = {}
    for i, (sid, v, b) in enumerate(scene.axons):
        if labeled_mask[i]:
            fp = footprint(v, scene.lo, na, z0, scene.crop_um, pix_um)
            fps[i] = fp
            total += fp
    return total, fps


def score_boutons(scene, labeled_mask, total, fps, na, z0,
                  axial_detect=AXIAL_DETECT, purity_thresh=PURITY_THRESH, pix_um=PIX_UM):
    """Return (detected, resolvable) bouton counts over labeled axons."""
    _, s_z = psf_sigmas_um(na)
    npx = total.shape[0]
    det = res = 0
    for i, (sid, v, b) in enumerate(scene.axons):
        if not labeled_mask[i] or len(b) == 0:
            continue
        aw = np.exp(-((b[:, 2] - z0) ** 2) / (2 * s_z ** 2))
        gx = ((b[:, 0] - scene.lo[0]) / pix_um).astype(int)
        gy = ((b[:, 1] - scene.lo[1]) / pix_um).astype(int)
        for j in range(len(b)):
            if aw[j] <= axial_detect or not (0 <= gx[j] < npx and 0 <= gy[j] < npx):
                continue
            det += 1
            tv = total[gx[j], gy[j]]
            pv = fps[i][gx[j], gy[j]]
            if tv > 0 and pv / tv >= purity_thresh:
                res += 1
    return det, res


def sweep(scene, densities, nas, rng, z0=None):
    """Sweep labeling density x NA. Returns a list of result dicts."""
    if z0 is None:
        z0 = scene.default_plane()
    rows = []
    for na in nas:
        for dens in densities:
            mask = rng.random(len(scene.axons)) < dens
            total, fps = image_labeled(scene, mask, na, z0)
            det, res = score_boutons(scene, mask, total, fps, na, z0)
            rows.append(dict(na=na, density=dens, labeled_axons=int(mask.sum()),
                             detected=det, resolvable=res,
                             resolvable_frac=(res / det if det else float("nan"))))
    return rows
