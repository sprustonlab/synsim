"""All tunable simulation knobs in one place, ready to drive a GUI / webpage.

`SimParams` is the single source of truth: every module sources its default
constants from `DEFAULTS`. Each field carries UI metadata (label, range, step,
unit, group) so a front-end can build sliders directly from `ui_spec()`, and
`to_dict`/`from_dict` give a clean JSON round-trip with the browser.
"""
from dataclasses import dataclass, field, asdict, fields


def _f(default, label, lo, hi, step, unit, group, desc):
    return field(default=default, metadata=dict(
        label=label, min=lo, max=hi, step=step, unit=unit, group=group, desc=desc))


@dataclass
class SimParams:
    # --- Microscope ---
    na: float = _f(0.6, "Numerical aperture (NA)", 0.3, 1.0, 0.05, "", "Microscope",
                   "Drives lateral (~1/NA) and axial (~1/NA^2) PSF; 0.6 = 2P-RAM mesoscope.")
    # --- Labeling ---
    density: float = _f(0.20, "Labeling density", 0.0, 1.0, 0.01, "fraction", "Labeling",
                        "Fraction of axons fluorescently labeled.")
    # --- Synapses ---
    bouton_spacing_um: float = _f(4.3, "Mean bouton spacing", 1.0, 10.0, 0.1, "um", "Synapses",
                                  "Sets synapse density; 4.3 um -> 0.23 syn/um (LICONN).")
    bouton_jitter: float = _f(0.2, "Spacing jitter", 0.0, 0.5, 0.05, "fraction", "Synapses",
                              "Random +/- fraction of spacing along the skeleton.")
    # --- Field of view (the imaged volume) ---
    fov_xy_um: float = _f(30.0, "FOV lateral (xy)", 5.0, 68.0, 1.0, "um", "Field of view",
                          "Lateral side of the imaged box.")
    fov_z_um: float = _f(19.0, "FOV axial (z)", 2.0, 19.0, 0.5, "um", "Field of view",
                         "Axial depth of the imaged box (LICONN volume is ~19 um deep).")
    edge_xy_um: float = _f(3.0, "Exclude edge (xy)", 0.0, 15.0, 0.5, "um", "Field of view",
                           "Don't score boutons within this xy margin (truncated PSF neighborhood).")
    edge_z_um: float = _f(5.0, "Exclude edge (z)", 0.0, 15.0, 0.5, "um", "Field of view",
                          "Don't score boutons within this z margin (axial PSF is large).")
    region: float = _f(0.0, "Region (sub-volume)", 0.0, 50.0, 1.0, "", "Field of view",
                       "Which neighborhood to image: 0 = volume center; other integers pick "
                       "a different reproducible random sub-volume.")
    # --- Imaging ---
    pix_um: float = _f(0.15, "Pixel size", 0.05, 0.5, 0.05, "um", "Imaging",
                       "Sim image pixel size (display sections only).")
    # --- Metric ---
    purity_thresh: float = _f(0.70, "Resolvability purity threshold", 0.5, 1.0, 0.05, "fraction",
                              "Metric",
                              "Parent axon's own bouton(s) must supply > this fraction of the bouton "
                              "signal at the focal plane (i.e. no other axon's bouton in the resolved spot).")

    def to_dict(self):
        return asdict(self)

    @classmethod
    def from_dict(cls, d):
        valid = {f.name for f in fields(cls)}
        return cls(**{k: v for k, v in d.items() if k in valid})


DEFAULTS = SimParams()


def ui_spec():
    """List of {name, value, label, min, max, step, unit, group, desc} for every
    knob - directly consumable by a front-end to build controls."""
    d = DEFAULTS
    return [dict(name=f.name, value=getattr(d, f.name), **f.metadata)
            for f in fields(SimParams)]
