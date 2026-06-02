"""Run the labeling-density x NA resolvability sweep and print a table.

    python scripts/run_sweep.py
"""
import numpy as np

from synsim import build_scene, sweep


def main():
    rng = np.random.default_rng(0)
    scene = build_scene(rng=rng)
    print(f"FOV {scene.fov[0]:.0f}x{scene.fov[1]:.0f}x{scene.fov[2]:.0f} um: "
          f"{len(scene.axons)} axons | length {scene.total_length_um:.0f} um | "
          f"boutons {scene.n_boutons} | synapses/um {scene.synapses_per_um:.3f}\n")

    rows = sweep(scene, densities=(0.05, 0.20, 0.60), nas=(0.4, 0.6, 0.8), rng=rng)
    print(f"{'NA':>4} {'dens':>5} {'lab.ax':>7} {'scored':>7} {'resolv':>7} {'resolv%':>8}")
    for r in rows:
        pct = 100 * r["resolvable_frac"]
        print(f"{r['na']:>4.1f} {r['density']:>5.2f} {r['labeled_axons']:>7} "
              f"{r['scored']:>7} {r['resolvable']:>7} {pct:>7.1f}%")


if __name__ == "__main__":
    main()
