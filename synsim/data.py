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


def voxel_downsample(v_um, ds_um):
    """Reduce a point set to one representative per ds_um voxel (uniform thinning)."""
    cells = np.unique(np.floor(v_um / ds_um).astype(np.int64), axis=0)
    return (cells + 0.5) * ds_um


def load_fov(kind="axons", fov_xy_um=30.0, fov_z_um=19.0, min_pts=50, data_dir=None):
    """Objects with >= min_pts vertices inside a 3D FOV box centered in the volume.

    Returns (objects, origin_xyz, fov_xyz) where objects is a list of
    (id, verts_um) holding only the in-box vertices, origin is the box min corner,
    and fov_xyz = (fov_xy, fov_xy, fov_z).
    """
    mn = np.array([np.inf] * 3)
    mx = np.array([-np.inf] * 3)
    sx = sy = n = 0.0
    for _, v in iter_objects(kind, data_dir):
        mn = np.minimum(mn, v.min(0))
        mx = np.maximum(mx, v.max(0))
        sx += v[:, 0].sum(); sy += v[:, 1].sum(); n += v.shape[0]
    cx, cy = sx / n, sy / n
    zc = 0.5 * (mn[2] + mx[2])
    fov = np.array([fov_xy_um, fov_xy_um, fov_z_um])
    origin = np.array([cx - fov_xy_um / 2, cy - fov_xy_um / 2, zc - fov_z_um / 2])
    hi = origin + fov
    out = []
    for sid, v in iter_objects(kind, data_dir):
        m = np.all((v >= origin) & (v < hi), axis=1)
        if m.sum() >= min_pts:
            out.append((sid, v[m]))
    return out, origin, fov
