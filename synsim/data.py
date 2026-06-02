"""Load extracted LICONN axon/dendrite geometry from the data/ cache.

The cache (data/axons, data/dends) is produced by scripts/extract_*.py: per-chunk
.npz with vertices in VOXEL coordinates (uint16) plus per-object offsets + ids.
Native voxel size is 18 x 18 x 24 nm.
"""
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
DATA_DIR = REPO / "data"
RES_NM = np.array([18.0, 18.0, 24.0])
UM = RES_NM / 1000.0                       # voxel -> um per axis


def iter_objects(kind="axons", data_dir=None):
    """Yield (id, verts_um Nx3) for every object in the cache."""
    base = Path(data_dir) if data_dir else DATA_DIR
    for f in sorted((base / kind).glob("chunk_*.npz")):
        d = np.load(f)
        verts, offs, ids = d["verts"], d["offsets"], d["ids"]
        for k, sid in enumerate(ids):
            yield int(sid), verts[offs[k]:offs[k + 1]].astype(np.float64) * UM


def load_crop(kind="axons", crop_um=20.0, min_pts=50, data_dir=None):
    """Objects with >= min_pts vertices inside a central lateral crop.

    Returns (axons, lo_xy, hi_xy) where axons is a list of (id, verts_um) holding
    only the in-crop vertices.
    """
    sx = sy = n = 0.0
    for _, v in iter_objects(kind, data_dir):
        sx += v[:, 0].sum(); sy += v[:, 1].sum(); n += v.shape[0]
    cx, cy = sx / n, sy / n
    lo = np.array([cx - crop_um / 2, cy - crop_um / 2])
    hi = lo + crop_um
    out = []
    for sid, v in iter_objects(kind, data_dir):
        m = (v[:, 0] >= lo[0]) & (v[:, 0] < hi[0]) & (v[:, 1] >= lo[1]) & (v[:, 1] < hi[1])
        if m.sum() >= min_pts:
            out.append((sid, v[m]))
    return out, lo, hi
