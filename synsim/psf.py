"""Anisotropic 2P-RAM point-spread function.

Anchored to Sofroniew et al. 2016 (eLife e14472): at NA 0.6 the measured PSF is
0.66 um lateral / 4.09 um axial FWHM. A single NA knob drives both axes,
lateral ~ 1/NA and axial ~ 1/NA^2.
"""
_FWHM_XY_NA06 = 0.66          # um, lateral FWHM at NA 0.6
_FWHM_Z_NA06 = 4.09           # um, axial FWHM at NA 0.6
_FWHM_TO_SIGMA = 1.0 / 2.3548200450309493
_REF_NA = 0.6


def psf_sigmas_um(na):
    """Return (sigma_xy, sigma_z) in um for the given numerical aperture."""
    fwhm_xy = _FWHM_XY_NA06 * (_REF_NA / na)
    fwhm_z = _FWHM_Z_NA06 * (_REF_NA / na) ** 2
    return fwhm_xy * _FWHM_TO_SIGMA, fwhm_z * _FWHM_TO_SIGMA
