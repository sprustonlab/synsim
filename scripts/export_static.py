"""Export the synsim dashboard as a standalone static site (no server, no Python).

Every parameter change in the live dashboard runs real compute, but the work
factors cleanly:

  - The optical-section image depends only on (na, density) at a fixed geometry.
  - A scored bouton's PURITY is independent of `purity_thresh` and of the edge
    margins (those only decide which boutons are counted). So if we store every
    imaged bouton's (x, y, z, purity), the browser can apply `edge_xy`, `edge_z`
    and `purity_thresh` live - and recompute the stats, histogram and bouton
    colors - with no server.

So we precompute a grid over (na x density) at the default geometry. Per cell we
write a normalized grayscale PNG of the central section plus a small JSON of its
imaged boutons. The static page (docs/index.html) looks up the nearest cell and
does the rest client-side.

    python scripts/export_static.py                 # -> docs/
    python scripts/export_static.py --out site      # custom output dir

Host the output dir on GitHub Pages (Settings -> Pages -> deploy from docs/).
"""
import argparse
import json
from pathlib import Path

import numpy as np
import matplotlib
from PIL import Image

from synsim import SimParams, ui_spec, build_scene, evaluate, section_image

# Knobs the viewer can change by indexing the precomputed grid.
NA_GRID = [round(0.30 + 0.05 * i, 2) for i in range(15)]          # 0.30 .. 1.00
DENSITY_GRID = [0.02, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30,
                0.35, 0.40, 0.45, 0.50, 0.60]

# Knobs the viewer can change for free, client-side (they only filter/threshold).
FREE_CONTROLS = ["edge_xy_um", "edge_z_um", "purity_thresh"]
GRID_CONTROLS = ["na", "density"]


def make_luts():
    luts = {}
    for name in ("gray", "viridis", "magma", "inferno", "hot", "cividis"):
        c = matplotlib.colormaps[name](np.linspace(0, 1, 256))[:, :3]
        luts[name] = (c * 255).astype(int).tolist()
    luts["green"] = [[0, i, 0] for i in range(256)]
    luts["red"] = [[i, 0, 0] for i in range(256)]
    return luts


def section_png(scene, mask, na, pix_um, path):
    """Write the central section as a normalized 8-bit grayscale PNG.

    Stored row-major over the (nx, ny) array, matching the live dashboard's flat
    index ix*ny+iy, so the browser reconstructs the exact same image.
    """
    img, _ = section_image(scene, mask, na, pix_um)                  # (nx, ny) float
    mx = float(img.max())
    norm = img / mx if mx > 0 else img
    u8 = np.clip(np.round(norm * 255), 0, 255).astype(np.uint8)      # (nx, ny)
    Image.fromarray(u8, mode="L").save(path, optimize=True)          # size = (ny, nx)
    return img.shape


def cell_boutons(scene, density, na):
    """All imaged boutons (edges off) with rel coords + purity, for one cell."""
    p = SimParams(na=na, density=density, edge_xy_um=0.0, edge_z_um=0.0)
    mask = np.random.default_rng(0).random(len(scene.axons)) < density
    r = evaluate(scene, mask, p)
    lo = scene.origin
    pos = r["positions"]
    if len(pos):
        bx = (pos[:, 0] - lo[0]).round(4).tolist()
        by = (pos[:, 1] - lo[1]).round(4).tolist()
        bz = (r["z_true"] - lo[2]).round(4).tolist()
        pur = r["purity"].round(4).tolist()
    else:
        bx = by = bz = pur = []
    return int(r["labeled_axons"]), {"bx": bx, "by": by, "bz": bz, "pur": pur}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="docs", help="output directory (default: docs)")
    args = ap.parse_args()

    out = Path(args.out)
    cells = out / "cells"
    cells.mkdir(parents=True, exist_ok=True)

    defaults = SimParams()
    print("building default scene ...", flush=True)
    scene = build_scene(rng=np.random.default_rng(0), params=defaults)
    nx = ny = None
    labeled_axons = [0] * len(DENSITY_GRID)

    total = len(NA_GRID) * len(DENSITY_GRID)
    done = 0
    for di, dens in enumerate(DENSITY_GRID):
        mask = np.random.default_rng(0).random(len(scene.axons)) < dens
        labeled_axons[di] = int(mask.sum())
        for ni, na in enumerate(NA_GRID):
            shape = section_png(scene, mask, na, defaults.pix_um,
                                cells / f"n{ni}_d{di}.png")
            nx, ny = int(shape[0]), int(shape[1])
            _, bdat = cell_boutons(scene, dens, na)
            (cells / f"n{ni}_d{di}.json").write_text(
                json.dumps(bdat), encoding="utf-8")
            done += 1
            if done % 10 == 0 or done == total:
                print(f"  {done}/{total} cells", flush=True)

    spec = {c["name"]: c for c in ui_spec()}
    meta = {
        "na_grid": NA_GRID,
        "density_grid": DENSITY_GRID,
        "nx": nx, "ny": ny,
        "fov": [float(scene.fov[0]), float(scene.fov[1])],
        "fov_z": float(scene.fov[2]),
        "pix_um": defaults.pix_um,
        "n_axons": len(scene.axons),
        "synapses_per_um": round(float(scene.synapses_per_um), 4),
        "labeled_axons": labeled_axons,
        "grid_controls": [spec[k] for k in GRID_CONTROLS],
        "free_controls": [spec[k] for k in FREE_CONTROLS],
        "fixed": {
            "fov_xy_um": defaults.fov_xy_um, "fov_z_um": defaults.fov_z_um,
            "bouton_spacing_um": defaults.bouton_spacing_um,
            "bouton_jitter": defaults.bouton_jitter, "region": defaults.region,
            "pix_um": defaults.pix_um,
        },
        "defaults": {k: getattr(defaults, k) for k in GRID_CONTROLS + FREE_CONTROLS},
    }
    (out / "meta.json").write_text(json.dumps(meta), encoding="utf-8")
    (out / "luts.json").write_text(json.dumps(make_luts()), encoding="utf-8")
    print(f"wrote {total} cells + meta to {out}/", flush=True)


if __name__ == "__main__":
    main()
