"""All tunable simulation knobs in one place, ready to drive a GUI / webpage.

`SimParams` is the single source of truth: every module sources its default
constants from `DEFAULTS`. Each field carries UI metadata (label, range, step,
unit, group) so a front-end can build sliders directly from `ui_spec()`, and
`to_dict`/`from_dict` give a clean JSON round-trip with the browser.
"""
from dataclasses import dataclass, field, asdict, fields
from typing import Optional


def _f(default, label, lo, hi, step, unit, group, desc, nullable=False):
    return field(default=default, metadata=dict(
        label=label, min=lo, max=hi, step=step, unit=unit, group=group,
        desc=desc, nullable=nullable))


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
                                  "Calibrated 4.3 um -> 0.23 synapses/um (LICONN).")
    bouton_jitter: float = _f(0.2, "Spacing jitter", 0.0, 0.5, 0.05, "fraction", "Synapses",
                              "Random +/- fraction of spacing along the skeleton.")
    # --- Volume / imaging ---
    crop_um: float = _f(20.0, "Crop size", 5.0, 60.0, 1.0, "um", "Volume",
                        "Lateral side of the simulated imaging field.")
    pix_um: float = _f(0.15, "Pixel size", 0.05, 0.5, 0.05, "um", "Imaging",
                       "Sim image pixel size.")
    z0_um: Optional[float] = _f(None, "Imaging plane z0", 0.0, 20.0, 0.5, "um", "Imaging",
                                "Depth of the imaged plane; null = densest bouton plane.",
                                nullable=True)
    # --- Metric thresholds ---
    axial_detect: float = _f(0.5, "Axial detection threshold", 0.1, 1.0, 0.05, "", "Metric",
                             "Bouton is in the optical section if axial PSF weight > this.")
    purity_thresh: float = _f(0.70, "Resolvability purity threshold", 0.5, 1.0, 0.05, "fraction",
                              "Metric",
                              "Parent axon must exceed this fraction of the pixel signal.")

    def to_dict(self):
        return asdict(self)

    @classmethod
    def from_dict(cls, d):
        valid = {f.name for f in fields(cls)}
        return cls(**{k: v for k, v in d.items() if k in valid})


DEFAULTS = SimParams()


def ui_spec():
    """List of {name, value, label, min, max, step, unit, group, desc, nullable}
    for every knob - directly consumable by a front-end to build controls."""
    d = DEFAULTS
    return [dict(name=f.name, value=getattr(d, f.name), **f.metadata)
            for f in fields(SimParams)]
