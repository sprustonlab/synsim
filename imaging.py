"""Shared imaging core: voxelize -> anisotropic 3D PSF convolution.

Single source of truth for the optical model, used by every figure/metric so a
sanity check on one validates all of them.
"""
import numpy as np
from scipy import ndimage

from spike_sim import psf_sigmas_um  # PSF sigma(NA), anchored to Sofroniew 2016


def voxelize(v_um, origin, shape, grid, dilate=1, fill=True):
    """Rasterize a point set (um) into a binary volume on a shared grid."""
    idx = np.floor((v_um - origin) / grid).astype(int)
    ok = np.all((idx >= 0) & (idx < np.array(shape)), axis=1)
    vol = np.zeros(shape, bool)
    idx = idx[ok]
    vol[idx[:, 0], idx[:, 1], idx[:, 2]] = True
    if dilate:
        vol = ndimage.binary_dilation(vol, iterations=dilate)
    if fill:
        vol = ndimage.binary_fill_holes(vol)
    return vol.astype(np.float32)


def blur3d(vol, na, grid):
    """Convolve a volume with the anisotropic 3D Gaussian PSF for a given NA."""
    s_xy, s_z = psf_sigmas_um(na)
    return ndimage.gaussian_filter(vol, (s_xy / grid, s_xy / grid, s_z / grid))


def section_index(z0_um, origin, grid):
    return int(round((z0_um - origin[2]) / grid))
