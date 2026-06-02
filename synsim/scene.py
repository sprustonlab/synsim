"""Assemble an imaging scene: cropped axons with boutons placed along each."""
import numpy as np

from .data import load_crop
from .boutons import skeleton_and_boutons


class Scene:
    """A cropped set of axons, each with vertices and placed boutons.

    axons: list of (id, verts_um Nx3, boutons_xyz Mx3).
    """

    def __init__(self, axons, lo, hi, total_length_um, crop_um):
        self.axons = axons
        self.lo = lo
        self.hi = hi
        self.total_length_um = total_length_um
        self.crop_um = crop_um

    @property
    def n_boutons(self):
        return sum(len(b) for _, _, b in self.axons)

    @property
    def synapses_per_um(self):
        return self.n_boutons / self.total_length_um if self.total_length_um else float("nan")

    def default_plane(self):
        """z (um) of the most populated bouton plane (+/-1 um neighborhood)."""
        bz = np.concatenate([b[:, 2] for _, _, b in self.axons if len(b)])
        return float(bz[np.argmax([(np.abs(bz - z) <= 1.0).sum() for z in bz])])


def build_scene(crop_um=20.0, rng=None, kind="axons", data_dir=None):
    """Load a crop and place boutons on every axon. Returns a Scene."""
    rng = rng if rng is not None else np.random.default_rng(0)
    axons, lo, hi = load_crop(kind, crop_um, data_dir=data_dir)
    out, total = [], 0.0
    for sid, v in axons:
        length, b = skeleton_and_boutons(v, rng)
        total += length
        out.append((sid, v, b))
    return Scene(out, lo, hi, total, crop_um)
