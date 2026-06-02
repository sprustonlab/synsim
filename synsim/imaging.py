"""Point-based optical-section imaging (the 'points approach').

A labeled axon's image at plane z0 is its vertices, each weighted by the axial
PSF around z0, accumulated laterally and blurred by the lateral PSF. This single
optical section was validated against a true 3D-convolution z0-slice (Pearson
r = 0.996), so it is used as the fast imaging operator.
"""
import numpy as np
from scipy import ndimage

from .psf import psf_sigmas_um
from .params import DEFAULTS

PIX_UM = DEFAULTS.pix_um      # sim image pixel size (default; override per call)


def footprint(v_um, lo_xy, na, z0, crop_um, pix_um=PIX_UM):
    """Axial-weighted, laterally-blurred image of one axon's points at plane z0."""
    s_xy, s_z = psf_sigmas_um(na)
    w = np.exp(-((v_um[:, 2] - z0) ** 2) / (2 * s_z ** 2))
    gx = ((v_um[:, 0] - lo_xy[0]) / pix_um).astype(int)
    gy = ((v_um[:, 1] - lo_xy[1]) / pix_um).astype(int)
    npx = int(round(crop_um / pix_um))
    ok = (gx >= 0) & (gx < npx) & (gy >= 0) & (gy < npx)
    img = np.zeros((npx, npx))
    np.add.at(img, (gx[ok], gy[ok]), w[ok])
    return ndimage.gaussian_filter(img, s_xy / pix_um)


def section_image(scene, labeled_mask, na, pix_um=PIX_UM):
    """Display image of the central optical section: all labeled axons' points,
    axial-weighted to the FOV-center plane, splatted in xy and lateral-blurred.

    Returns (img [nx, ny], z_plane). One pass + one blur (fast for the dashboard).
    """
    s_xy, s_z = psf_sigmas_um(na)
    lo = scene.origin
    fov = scene.fov
    z_plane = lo[2] + fov[2] / 2.0
    nx = int(round(fov[0] / pix_um))
    ny = int(round(fov[1] / pix_um))
    img = np.zeros((nx, ny))
    for i, (sid, v, b) in enumerate(scene.axons):
        if not labeled_mask[i]:
            continue
        w = np.exp(-((v[:, 2] - z_plane) ** 2) / (2 * s_z ** 2))
        gx = ((v[:, 0] - lo[0]) / pix_um).astype(int)
        gy = ((v[:, 1] - lo[1]) / pix_um).astype(int)
        ok = (gx >= 0) & (gx < nx) & (gy >= 0) & (gy < ny)
        np.add.at(img, (gx[ok], gy[ok]), w[ok])
    return ndimage.gaussian_filter(img, s_xy / pix_um), z_plane
