"""synsim - simulation of sparse 2P axon imaging against LICONN ground truth.

Quantifies how often in-vivo 2P imaging of sparsely labeled axons can resolve
individual axons (and thus assign their synapses), as a function of labeling
density and numerical aperture, using the real LICONN axon population as truth.
"""
from .params import SimParams, DEFAULTS, ui_spec
from .psf import psf_sigmas_um
from .data import load_fov, iter_objects, voxel_downsample, RES_NM, UM
from .boutons import skeleton_and_boutons, skeleton_points, BOUTON_SPACING_UM
from .imaging import footprint, section_image, PIX_UM
from .scene import Scene, build_scene
from .metrics import evaluate, operating_point, sweep, PURITY_THRESH

__version__ = "0.2.0"
__all__ = [
    "SimParams", "DEFAULTS", "ui_spec",
    "psf_sigmas_um", "load_fov", "iter_objects", "voxel_downsample", "RES_NM", "UM",
    "skeleton_and_boutons", "skeleton_points", "BOUTON_SPACING_UM",
    "footprint", "section_image", "PIX_UM", "Scene", "build_scene",
    "evaluate", "operating_point", "sweep", "PURITY_THRESH",
]
