"""Run the labeling-density x NA resolvability sweep and print a table.

    python scripts/run_sweep.py
"""
import numpy as np

from synsim import build_scene, sweep


def main():
    rng = np.random.default_rng(0)
    scene = build_scene(crop_um=20.0, rng=rng)
    print(f"crop {scene.crop_um:.0f} um: {len(scene.axons)} axons | "
          f"length {scene.total_length_um:.0f} um | boutons {scene.n_boutons} | "
          f"synapses/um {scene.synapses_per_um:.3f} (target 0.233)")
    z0 = scene.default_plane()
    print(f"imaging plane z0 = {z0:.1f} um\n")

    rows = sweep(scene, densities=(0.05, 0.20, 0.60), nas=(0.4, 0.6, 0.8),
                 rng=rng, z0=z0)
    print(f"{'NA':>4} {'dens':>5} {'lab.ax':>7} {'det':>6} {'resolv':>7} {'resolv%':>8}")
    for r in rows:
        pct = 100 * r["resolvable_frac"]
        print(f"{r['na']:>4.1f} {r['density']:>5.2f} {r['labeled_axons']:>7} "
              f"{r['detected']:>6} {r['resolvable']:>7} {pct:>7.1f}%")


if __name__ == "__main__":
    main()
