# Sparse axon-imaging simulation — plan

Status: **core built and validated; this doc now specifies the GUI.**
The imaging model, bouton placement, and metrics live in the `synsim` package and
are validated on real LICONN geometry. What remains is the front-end (controls +
panels) and the data-baking step that feeds it.

---

## A. What the simulation argues

**Scientific question (from the Magee/Spruston thread):** can in-vivo 2P imaging of
*sparsely labeled* CA3 axons resolve individual axons well enough to do the post-hoc
connection step — or does **axon bundling / co-fasciculation** (Jeff's stated failure
mode, believed rare but with *no data*) merge neighbors into unresolvable blobs? And
what is the **yield vs. purity** trade-off: sparser labeling buys separability but costs
candidate synapses.

**Why it's credible:** the **real LICONN axon population** (18,667 axons, `data/axons/`)
is ground truth; we impose a realistic 2P-RAM imaging model on top. Jeff's hand-wave
becomes a number.

---

## B. Imaging model (implemented in `synsim`)

- **Ground truth:** real axon geometry, LICONN volume 68.8 x 87.7 x 19.1 um, native
  18 x 18 x 24 nm voxels. Dendrites (1,643) also cached for a future postsynaptic view.
- **PSF (`synsim.psf`):** anisotropic Gaussian, anchored to Sofroniew 2016 (eLife e14472):
  NA 0.6 -> lateral FWHM 0.66 um, axial FWHM 4.09 um. One NA knob drives both:
  lateral ~ 1/NA, axial ~ 1/NA^2.
- **Sparse labeling:** label a random fraction `density` of axons; only labeled axons
  fluoresce (mirrors viral-titer / birthdate-matched sparsity).
- **Bouton placement (`synsim.boutons`):** each axon is 3D-skeletonized; boutons are
  dropped by **walking the skeleton** (MST + DFS, branch-aware) every `bouton_spacing_um`
  +/- `bouton_jitter`. Default spacing 4.3 um reproduces the LICONN density
  (0.95 pre/um^3 / 4.084 um-length/um^3 = 0.233 syn/um); spacing is a GUI knob.
- **FOV + edge exclusion:** imaging is restricted to a box `fov_xy_um` x `fov_xy_um`
  x `fov_z_um`. Boutons within `edge_xy_um` / `edge_z_um` of the boundary are NOT scored
  (their PSF neighborhood would be truncated by the FOV); axons in that margin still
  contribute signal, so interior boutons see their true neighborhood.
- **Imaging one central plane (`synsim.metrics`):** a single optical section at the FOV
  center in z. Each interior bouton is read at that plane: the signal is evaluated at
  (x_bouton, y_bouton, z_center) - lateral at the bouton, axial at the focal plane - so a
  bouton off the central plane is defocused, as in a real single-plane acquisition.
- **Metric (gridless):** a bouton is **resolvable** if its parent axon supplies
  > `purity_thresh` of the PSF-weighted signal at its plane location. Computed by scaling
  coordinates so the PSF is isotropic, KD-tree over labeled axon points, summing
  exp(-d^2/2) over neighbors within 3 sigma (parent vs. total). Every interior bouton is
  scored; yield = labeled axons; the trade-off curve is resolvable-vs-density.

## Synapse ground truth (we generate it; not in the public data)

The public LICONN release has no synapse/partner annotations. We place boutons along
real axons (above), and reduce assignability to **"a synapse is assignable iff its parent
axon is resolvable at that bouton."** Boutons give real 3D synapse positions to score and
count; the overlap question reduces to axon separability.

---

## C. Current state (done)

- **Extraction:** all 18,667 axons + 1,643 dendrites pulled from the public GCS bucket
  to `data/` (gitignored). Scripts: `scripts/extract_axons.py`, `extract_dends.py`.
- **Package `synsim`:** `psf, data, boutons, imaging, scene, metrics, params`
  (editable-installed in `.venv-liconn`).
- **Validation:** along-axon density calibrated; resolvability monotonic in density and
  NA; 2D point model == 3D conv (r=0.996); arc-length placement gives quasi-regular
  boutons.
- **Params exposed:** every knob is a `SimParams` field with UI metadata;
  `ui_spec()` + `to_dict/from_dict` + `scripts/dump_params.py` give the JSON contract for
  the front-end; `metrics.operating_point(scene, params, rng)` = one params in, one
  result dict out.

---

## D. The GUI

### Controls (each maps 1:1 to a `SimParams` field; ranges from `ui_spec()`)

| Group | Control | Range (default) | Drives |
|-------|---------|-----------------|--------|
| Microscope | **NA** | 0.30–1.00 (0.60) | PSF (both axes) → panels 1,2,3,4 |
| Labeling | **Labeling density** | 0–1 (0.20) | which axons fluoresce → 1,3,4 |
| Synapses | **Bouton spacing (um)** | 1.0–10.0 (4.3) | synapse density → 1,3,4 |
| Synapses | **Spacing jitter** | 0–0.5 (0.20) | placement regularity → 1 |
| Field of view | **FOV lateral xy (um)** | 5–68 (30) | imaged box (xy) → 1,3 |
| Field of view | **FOV axial z (um)** | 2–19 (19) | imaged box depth → 1,2,3 |
| Field of view | **Exclude edge xy (um)** | 0–15 (3) | drop boundary boutons → 3,4 |
| Field of view | **Exclude edge z (um)** | 0–15 (5) | drop boundary boutons → 3,4 |
| Field of view | **Region (sub-volume)** | 0–50 (0) | 0=center, else random neighborhood → 1,3,4 |
| Imaging | **Pixel size (um)** | 0.05–0.5 (0.15) | display section sampling → 1 |
| Metric | **Purity threshold** | 0.5–1.0 (0.70) | resolvable cutoff → 1,3,4 |

Imaging is a **single optical section at the FOV center in z**; boutons off that plane are
defocused. A small **"Experiment context"** group will be added for the readout scaling
(P_connect 0.10/0.20/0.25; cells/cohort) — annotation only, not physics.

### Panels

1. **Microscope view (hero) — z0 optical section.** Live render of the labeled axons'
   imaged plane at current settings; each in-plane bouton ringed **green (resolvable)** or
   **red (confused)**; toggle ground-truth overlay. This is the interactive version of the
   vignette figure.
2. **Side view (xz) — depth of field.** Orthogonal projection showing the axial PSF smear
   at current NA / z0 — *why* axons at different depths merge into one plane. Makes the
   1/NA^2 axial term tangible; updates with NA and z0.
3. **Yield vs. purity curve.** Sweep labeling density (x-axis) at the current NA/spacing;
   plot resolvable-axon yield and resolvable fraction; mark the current operating point;
   annotate reference markers (Jeff ~250 synapses/cell; proposal 10–200/cell; cohort
   1,100 / 4,400 / 25,000).
4. **Readout / scoreboard.** Live numbers from `operating_point`: labeled axons, detected
   boutons, resolvable count + %, synapses/um; plus an experiment-scaled estimate of
   *usable* classified synapses (× P_connect × cells) shown next to Jeff's 250/cell sketch.

---

## E. Architecture & data baking

**Hybrid (recommended).** Ship a **decimated subvolume** (axon centerlines / sparse point
clouds, target < 5 MB) so the browser can recompute the microscope view live in JS/canvas
(splat axial-weighted points → separable Gaussian blur → re-score boutons). **Precompute**
the population yield/purity curves (density sweep at a few NAs) in Python → JSON, with the
live operating point overlaid. Honest about which numbers are live vs. precomputed.

- **Baking:** reduce the 192 MB vertex cache to per-axon centerlines (reuse
  `synsim.boutons.skeleton_points`) for one representative subvolume; emit compact JSON.
  Validate baked metrics vs. full-resolution Python (sanity gate).
- **Fallback if JS is too slow:** precompute-only — Python bakes a grid of rendered
  sections + results; sliders snap to the grid. (Spike showed ~1 s for 2,169 axons in
  Python, so live JS on a subvolume should be fine.)
- **Backend option (alternative):** a thin FastAPI server calling `operating_point` for
  fully continuous, exact sliders — drop the baking, lose static shareability.

---

## F. Milestones

- **M0 — core + validation.** DONE (`synsim`, spike, r=0.996, arc-length, params).
- **M1 — bake subvolume asset** (centerlines → JSON) + validate baked vs full-res.
- **M2 — static HTML shell** with controls from `dump_params.py` + Panel 1 (microscope
  view) rendering live in JS.
- **M3 — Panels 2–4** (side view, yield/purity curve, readout) + reference annotations.
- **M4 — polish**; optional postsynaptic extension using `data/dends/`.

---

## G. Open decisions
1. **Architecture:** hybrid static page (recommended) vs. thin Python backend.
2. **Representative subvolume** for the live view: the central crop (current default) or a
   specific neighborhood (e.g. a dense bundle to showcase the failure mode).
3. **Experiment-context knobs** in the readout (P_connect, cells/cohort) — include now or
   later?
