"""Place presynaptic boutons along an axon centerline by arc length.

Linear density is calibrated to the LICONN paper: 0.95 pre-synapses/um^3 divided
by 4.084 um axon-length/um^3 = 0.233 synapses/um, i.e. mean bouton spacing ~4.3 um
(en-passant range). Density therefore falls out of real axon length with no
dependence on the dendrite side.

Each axon is voxelized and 3D-skeletonized to recover its centerline. The
skeleton (which may branch) is turned into a tree (MST over neighbouring voxels);
boutons are dropped by walking that tree and placing one every `spacing_um`
(+/- jitter) of accumulated arc length, so spacing is quasi-regular like real
en-passant boutons rather than randomly clumped.
"""
import numpy as np
from scipy import ndimage
from scipy.spatial import cKDTree
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import minimum_spanning_tree, connected_components
from skimage.morphology import skeletonize

from .params import DEFAULTS

SKEL_UM = 0.25                            # voxel size for skeletonization
BOUTON_SPACING_UM = DEFAULTS.bouton_spacing_um   # 0.95 / 4.084 = 0.233 synapses/um
JITTER = DEFAULTS.bouton_jitter                  # +/- fraction of spacing


def skeleton_points(v_um, skel_um=SKEL_UM):
    """Voxelize an axon's points and return its 3D skeleton centerline (um)."""
    mn = v_um.min(0)
    idx = np.floor((v_um - mn) / skel_um).astype(int)
    shp = idx.max(0) + 3
    vol = np.zeros(shp, bool)
    vol[idx[:, 0] + 1, idx[:, 1] + 1, idx[:, 2] + 1] = True
    vol = ndimage.binary_dilation(vol, iterations=1)
    vol = ndimage.binary_fill_holes(vol)
    return (np.argwhere(skeletonize(vol)) - 1) * skel_um + mn


def _walk_and_place(skel, rng, spacing, skel_um, jitter):
    """Drop boutons every `spacing` (+/- jitter) of arc length along the skeleton tree."""
    n = len(skel)
    if n < 2:
        return skel.copy()
    pairs = cKDTree(skel).query_pairs(r=skel_um * 1.8, output_type="ndarray")
    if len(pairs) == 0:
        return skel[[rng.integers(n)]]
    w = np.linalg.norm(skel[pairs[:, 0]] - skel[pairs[:, 1]], axis=1)
    rows = np.concatenate([pairs[:, 0], pairs[:, 1]])
    cols = np.concatenate([pairs[:, 1], pairs[:, 0]])
    g = csr_matrix((np.concatenate([w, w]), (rows, cols)), shape=(n, n))
    mst = minimum_spanning_tree(g)
    mst = mst + mst.T
    coo = mst.tocoo()
    adj = [[] for _ in range(n)]
    for a, b, wt in zip(coo.row, coo.col, coo.data):
        adj[a].append((b, wt))
    deg = np.array([len(a) for a in adj])
    _, labels = connected_components(g, directed=False)

    def target():
        return spacing * (1 + rng.uniform(-jitter, jitter))

    chosen = []
    for comp in range(labels.max() + 1):
        nodes = np.where(labels == comp)[0]
        if len(nodes) < 2:
            continue
        leaves = nodes[deg[nodes] == 1]
        start = int(rng.choice(leaves)) if len(leaves) else int(nodes[0])
        stack = [(start, -1, 0.0, target())]      # (node, parent, accum, target)
        while stack:
            u, parent, accum, tgt = stack.pop()
            for v, wt in adj[u]:
                if v == parent:
                    continue
                a, t = accum + wt, tgt
                if a >= t:
                    chosen.append(v)
                    a -= t
                    t = target()
                stack.append((v, u, a, t))
    if not chosen:
        return skel[[rng.integers(n)]]
    return skel[np.array(chosen)]


def skeleton_and_boutons(v_um, rng, spacing_um=BOUTON_SPACING_UM, skel_um=SKEL_UM,
                         jitter=JITTER):
    """Return (length_um, boutons_xyz Nx3) for one axon's points."""
    skel = skeleton_points(v_um, skel_um)
    if len(skel) == 0:
        return 0.0, np.zeros((0, 3))
    length = len(skel) * skel_um
    return length, _walk_and_place(skel, rng, spacing_um, skel_um, jitter)
