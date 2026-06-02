"""Place presynaptic boutons along an axon centerline.

Linear density is calibrated to the LICONN paper: 0.95 pre-synapses/um^3 divided
by 4.084 um axon-length/um^3 = 0.233 synapses/um, i.e. mean bouton spacing ~4.3 um
(en-passant range). Density therefore falls out of real axon length with no
dependence on the dendrite side.

Each axon is voxelized and 3D-skeletonized to recover its centerline + length;
boutons are sampled along the skeleton at the calibrated spacing.
"""
import numpy as np
from scipy import ndimage
from skimage.morphology import skeletonize

SKEL_UM = 0.25               # voxel size for skeletonization
BOUTON_SPACING_UM = 4.3      # 0.95 / 4.084 = 0.233 synapses/um


def skeleton_and_boutons(v_um, rng, spacing_um=BOUTON_SPACING_UM, skel_um=SKEL_UM):
    """Return (length_um, boutons_xyz Nx3) for one axon's points."""
    mn = v_um.min(0)
    idx = np.floor((v_um - mn) / skel_um).astype(int)
    shp = idx.max(0) + 3
    vol = np.zeros(shp, bool)
    vol[idx[:, 0] + 1, idx[:, 1] + 1, idx[:, 2] + 1] = True
    vol = ndimage.binary_dilation(vol, iterations=1)
    vol = ndimage.binary_fill_holes(vol)
    skel = skeletonize(vol)
    sp = np.argwhere(skel)
    if sp.shape[0] == 0:
        return 0.0, np.zeros((0, 3))
    length = sp.shape[0] * skel_um
    skel_pts = (sp - 1) * skel_um + mn
    nb = max(1, int(round(length / spacing_um)))
    pick = rng.choice(sp.shape[0], size=min(nb, sp.shape[0]), replace=False)
    return length, skel_pts[pick]
