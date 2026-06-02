"""Resolvability metric: gridless, per-bouton, in 3D.

The axons flash (GCaMP), so what matters for assigning a synapse to a
presynaptic axon is whether MORE THAN ONE BOUTON sits in the same resolved spot -
not how much axon membrane is nearby. The emitters are therefore the BOUTONS
(each tagged by its parent axon), not the axon membrane.

We image a SINGLE optical section at the center of the FOV in z. Every labeled
bouton inside the FOV interior (FOV minus edge margins) is read at that plane:
the bouton-signal is evaluated at (x_bouton, y_bouton, z_center) - lateral at the
bouton, axial at the focal plane (so a bouton off the central plane is defocused).
A bouton is RESOLVABLE if its parent axon's own bouton(s) supply more than
`purity_thresh` of the PSF-weighted bouton-signal there, i.e. no other axon's
bouton falls within the resolved spot.

Implementation: scale coordinates by (1/sigma_xy, 1/sigma_xy, 1/sigma_z) so the
PSF is isotropic, build a KD-tree over labeled BOUTONS, and for each scored
bouton sum exp(-d^2/2) over neighbors within 3 sigma, splitting same-axon (parent)
vs. all boutons.
"""
import numpy as np
from scipy.spatial import cKDTree

from .psf import psf_sigmas_um
from .params import DEFAULTS

PURITY_THRESH = DEFAULTS.purity_thresh
_BALL = 3.0                    # neighbor cutoff in units of sigma (scaled space)


def evaluate(scene, labeled_mask, params):
    """Score all interior boutons of labeled axons. Returns a result dict with
    counts and per-bouton arrays (positions, purity, resolvable)."""
    s_xy, s_z = psf_sigmas_um(params.na)
    scale = np.array([1.0 / s_xy, 1.0 / s_xy, 1.0 / s_z])

    # emitter point cloud = all labeled BOUTONS (scaled), tagged by parent axon
    pts, tags = [], []
    for i, (sid, v, b) in enumerate(scene.axons):
        if labeled_mask[i] and len(b):
            pts.append(b * scale)
            tags.append(np.full(len(b), i))
    empty = dict(labeled_axons=int(np.sum(labeled_mask)), scored=0, resolvable=0,
                 resolvable_frac=float("nan"), positions=np.zeros((0, 3)),
                 purity=np.zeros(0), resolved=np.zeros(0, bool))
    if not pts:
        return empty
    V = np.concatenate(pts)
    A = np.concatenate(tags)
    tree = cKDTree(V)

    # single imaging plane (optical section) at the FOV center in z
    lo = scene.origin
    fov = scene.fov
    z_plane = lo[2] + fov[2] / 2.0

    # interior boutons of labeled axons, all read at the central plane
    ex, ez = params.edge_xy_um, params.edge_z_um
    bpos, bpar = [], []
    for i, (sid, v, b) in enumerate(scene.axons):
        if not labeled_mask[i] or len(b) == 0:
            continue
        rel = b - lo
        inside = ((rel[:, 0] >= ex) & (rel[:, 0] <= fov[0] - ex) &
                  (rel[:, 1] >= ex) & (rel[:, 1] <= fov[1] - ex) &
                  (rel[:, 2] >= ez) & (rel[:, 2] <= fov[2] - ez))
        for bb in b[inside]:
            bpos.append(np.array([bb[0], bb[1], z_plane]) * scale)   # read at the focal plane
            bpar.append(i)
    if not bpos:
        return empty
    bpos = np.array(bpos)
    bpar = np.array(bpar)

    neigh = tree.query_ball_point(bpos, r=_BALL)
    purity = np.zeros(len(bpos))
    for k, nb in enumerate(neigh):
        if not nb:
            continue
        nb = np.asarray(nb)
        w = np.exp(-0.5 * ((V[nb] - bpos[k]) ** 2).sum(1))
        tot = w.sum()
        if tot > 0:
            purity[k] = w[A[nb] == bpar[k]].sum() / tot
    resolved = purity >= params.purity_thresh
    scored = len(bpos)
    return dict(labeled_axons=int(np.sum(labeled_mask)), scored=scored,
                resolvable=int(resolved.sum()),
                resolvable_frac=float(resolved.mean()),
                positions=bpos / scale, purity=purity, resolved=resolved)


def operating_point(scene, params, rng):
    """Image + score one (na, density) operating point. Returns a GUI-ready dict."""
    mask = rng.random(len(scene.axons)) < params.density
    r = evaluate(scene, mask, params)
    return dict(na=params.na, density=params.density,
                labeled_axons=r["labeled_axons"], scored=r["scored"],
                resolvable=r["resolvable"], resolvable_frac=r["resolvable_frac"])


def sweep(scene, densities, nas, rng, params=DEFAULTS):
    """Sweep labeling density x NA. Returns a list of result dicts."""
    rows = []
    for na in nas:
        p = SimParams_with(params, na=na)
        for dens in densities:
            mask = rng.random(len(scene.axons)) < dens
            r = evaluate(scene, mask, p)
            rows.append(dict(na=na, density=dens, labeled_axons=r["labeled_axons"],
                             scored=r["scored"], resolvable=r["resolvable"],
                             resolvable_frac=r["resolvable_frac"]))
    return rows


def SimParams_with(params, **changes):
    """Return a copy of params with fields overridden."""
    d = params.to_dict()
    d.update(changes)
    return type(params).from_dict(d)
