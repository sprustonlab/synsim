"""Assemble an imaging scene: axons in a 3D FOV box, each with placed boutons.

For each axon we keep:
  - verts_ds: voxel-downsampled vertices (the fluorescent point cloud used to
    compute the 3D-PSF-weighted signal; uniform thinning leaves purity ratios
    unbiased while keeping the neighbor search fast).
  - boutons:  presynaptic sites placed along the axon (the things we score).
"""
import numpy as np

from .data import load_fov, voxel_downsample
from .boutons import skeleton_and_boutons
from .params import DEFAULTS

DS_UM = 0.25                  # downsample spacing for the signal point cloud


class Scene:
    def __init__(self, axons, origin, fov, total_length_um):
        self.axons = axons                  # list of (id, verts_ds Nx3, boutons Mx3)
        self.origin = np.asarray(origin)     # box min corner (xyz, um)
        self.fov = np.asarray(fov)           # (fov_xy, fov_xy, fov_z) um
        self.total_length_um = total_length_um

    @property
    def n_boutons(self):
        return sum(len(b) for _, _, b in self.axons)

    @property
    def synapses_per_um(self):
        return self.n_boutons / self.total_length_um if self.total_length_um else float("nan")


def build_scene(rng=None, params=None, kind="axons", data_dir=None, ds_um=DS_UM):
    """Load a FOV box and place boutons on every axon, per `params`. Returns a Scene."""
    params = params if params is not None else DEFAULTS
    rng = rng if rng is not None else np.random.default_rng(0)
    objs, origin, fov = load_fov(kind, params.fov_xy_um, params.fov_z_um, data_dir=data_dir)
    out, total = [], 0.0
    for sid, v in objs:
        length, b = skeleton_and_boutons(v, rng,
                                         spacing_um=params.bouton_spacing_um,
                                         jitter=params.bouton_jitter)
        total += length
        out.append((sid, voxel_downsample(v, ds_um), b))
    return Scene(out, origin, fov, total)
