"""synsim - simulation of sparse 2P axon imaging against LICONN ground truth.

Quantifies how often in-vivo 2P imaging of sparsely labeled axons can resolve
individual axons (and thus assign their synapses), as a function of labeling
density and numerical aperture, using the real LICONN axon population as truth.
"""
from .psf import psf_sigmas_um
from .data import load_crop, iter_objects, RES_NM, UM
from .boutons import skeleton_and_boutons, BOUTON_SPACING_UM
from .imaging import footprint, PIX_UM
from .scene import Scene, build_scene
from .metrics import (image_labeled, score_boutons, sweep,
                      AXIAL_DETECT, PURITY_THRESH)

__version__ = "0.1.0"
__all__ = [
    "psf_sigmas_um", "load_crop", "iter_objects", "RES_NM", "UM",
    "skeleton_and_boutons", "BOUTON_SPACING_UM", "footprint", "PIX_UM",
    "Scene", "build_scene", "image_labeled", "score_boutons", "sweep",
    "AXIAL_DETECT", "PURITY_THRESH",
]
