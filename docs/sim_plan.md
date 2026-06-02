# Sparse axon-imaging simulation — plan

Status: **proposal for review.** Nothing here is built yet. The spike (Part B) is
designed but not run; approval gates both the spike and the full build (Part C).

---

## A. What the simulation argues, and why

**Scientific question (from the Magee/Spruston email thread):** can in-vivo 2P imaging
of *sparsely labeled* CA3 axons actually resolve individual axons well enough to do the
post-hoc connection step — or does **axon bundling / co-fasciculation** (Jeff's stated
failure mode, which he believes is rare but has *no data* on) merge neighbors into
unresolvable blobs? And what is the **yield vs. purity** trade-off: sparser labeling
buys separability but costs candidate synapses.

**Why this is credible here:** we use the **real LICONN axon population** (18,667 axons,
already extracted to `data/axons/`) as ground truth, and impose a realistic 2P-RAM
imaging model on top. This turns Jeff's hand-wave ("bundling is rare") into a number.

**Two deliverables the sim produces:**
1. Axon-confusion vs. labeling density (the bundling question).
2. Yield vs. purity curve, annotated against the proposal's own yield numbers
   (Jeff's ~250 synapses/cell sketch; the proposal's 10–200/cell; cohort 1,100–25,000).

---

## Imaging model (the physics the sim implements)

- **Ground truth:** real axon geometry in the LICONN volume
  68.8 x 87.7 x 19.1 um, native 18 x 18 x 24 nm voxels.
- **PSF:** anisotropic 3D Gaussian. Anchored to Sofroniew 2016 (eLife e14472):
  at **NA 0.6**, lateral FWHM = **0.66 um**, axial FWHM = **4.09 um**.
- **NA slider (single knob, drives both axes):**
  - lateral FWHM(NA) = 0.66 * (0.6 / NA)        (lateral ~ 1/NA)
  - axial   FWHM(NA) = 4.09 * (0.6 / NA)^2       (axial ~ 1/NA^2)
  - sigma = FWHM / 2.3548. At NA 0.6: sigma_xy = 0.28 um, sigma_z = 1.74 um.
- **Sparse labeling:** label a random fraction `p` of axons (the "labeling density"
  knob). Only labeled axons fluoresce — mirrors viral-titer / birthdate-matched sparsity.
- **Imaging a plane (2-plane mesoscope):** the detected image at depth z0 is the
  labeled signal **axially weighted by the PSF around z0** (so axons within ~sigma_z of
  the plane pile in — this is where the NA knob bites), summed in z, then **laterally
  blurred** by the in-plane PSF. Depth/plane is a knob; optional edge-of-FOV mode swaps
  to the degraded PSF (lateral 0.89 / axial 6.88 um).

## Synapse ground truth (NOT in the public data — we generate it)

The public LICONN release has **no synapse / spine / partner annotations** — only image,
dense segmentation, and the axon/dendrite id lists. The paper's 71,269 spines and
pre->post assignments were not released.

**Chosen approach: boutons placed ALONG axons (statistical, length-driven).** Follow each
real axon's centerline and place presynaptic boutons as a point process whose linear
density reproduces the published volumetric synapse density. A synapse is then
**assignable iff its parent axon is resolvable at that bouton's location** (the axon-only
reduction). The overlap question reduces to axon separability, but boutons give us real
3D synapse positions to score and to count yield.

**Calibration (from the LICONN paper):**
- pre-synapse density 0.95 /um^3; axon-length density 4.084 um/um^3 (342.3 mm / 83,825 um^3).
- => synapses per um of axon = 0.95 / 4.084 = **0.233 /um  (mean bouton spacing ~4.3 um)**.
- Default: place boutons along each axon centerline at mean spacing 4.3 um (Poisson, or
  uniform + jitter); spacing is a knob. This yields the right volumetric density *given
  real axon length*, with no dependence on the dendrite side.

**Procedure per axon:** voxelize vertices -> 3D skeletonize -> centerline + length ->
place boutons at calibrated spacing along the centerline.

The spike (Part B) validates this: that skeleton lengths + bouton placement reproduce the
target ~0.23 synapses/um and a sane total synapse count in the crop.

## Confusion / yield / purity metrics (candidate definitions)

- **Confusion (impurity):** fraction of detected fluorescence mass landing in **mixed
  pixels** — pixels where no single labeled axon contributes > X% of the signal
  (default X = 70%). Rises with density and as NA falls.
- **Resolvable fraction (purity):** fraction of labeled axons that own a contiguous set
  of pixels they dominate (i.e., separable from neighbors). 1 - (merged fraction).
- **Yield:** number of labeled axons present in the imaged volume, scaled to "candidate
  presynaptic partners" and then to candidate synapses via P_connect — directly
  comparable to Jeff's ~250/cell.
- Recommended headline metric: **resolvable-axon yield = (labeled axons) x (purity)** —
  captures both halves of the trade-off in one curve.

---

## B. The SPIKE (proof-of-concept, before any full build)

**Goal:** prove the whole core works on real geometry and is sensible + tractable —
*and* tell us which build architecture (Part C) is viable.

**What it does (one Python script, isolated uv env, target < 2 min):**
1. Load axons intersecting a small lateral crop (~30 x 30 um), one z-plane.
2. **Skeleton + bouton check:** voxelize -> 3D skeletonize each axon -> length; place
   boutons at mean spacing 4.3 um; verify synapses/um ~= 0.23 and total count is sane.
3. Voxelize axons onto a sim grid (150 nm pixels), apply the anisotropic PSF at NA 0.6.
4. **Per-bouton resolvability:** at each bouton's pixel, parent-axon signal fraction
   (purity); assignable if purity > 70%. Sweep labeling density (5/20/60%) x NA (0.4/0.6/0.8).
5. Emit: (a) table of yield / resolvable-fraction / synapses-per-um, (b) one PNG showing
   two real axons merging as density/blur rises, with boutons marked.

**Decision criteria the spike must satisfy before we build:**
- Metrics are **monotonic and physically sensible** (confusion up with density; up as
  NA down).
- The PSF/voxelize/metric pipeline runs **fast enough** to judge browser feasibility
  (informs A vs B vs C below).
- The demo PNG **visibly shows** real-axon merging — i.e., the bundling effect is real
  in this dataset, not an artifact of my model.

**Cost if approved:** ~30 min of my time; one throwaway script + one PNG. No commitment
to the full build.

---

## B-results. Spike outcome (PASSED — `spike_sim.py`, `spike_demo.png`)

- Crop 20x20 um: 2,169 real axons; total length 19,147 um; 4,742 boutons ->
  **0.248 synapses/um** (target 0.233) — along-axon bouton model is calibrated.
- Resolvable% **monotonic** both ways: NA 0.6 -> 35/2/0% at 5/20/60% labeling;
  at 5% labeling -> 12/35/78% for NA 0.4/0.6/0.8.
- Yield-vs-purity tension is real: lower NA detects MORE boutons/plane (155 vs 46)
  but resolves fewer. Headline: at NA 0.6, ~35% resolvable at 5% labeling, ~2% at 20%.
- Runtime 11.5 s for 2,169 axons x 9 conditions -> **browser-live feasible** for one
  subvolume; **architecture C (hybrid) confirmed viable.**
- Caveats (build-time knobs, not ground truth): purity threshold 0.70, pixel 150 nm,
  axial-detect cutoff 0.5; boutons uniform along skeleton; single plane / crop / seed.

## C. The FULL BUILD (proposed, gated on a good spike)

### Architecture options
- **(A) Precompute-only:** Python sweeps the NA x density (x depth) grid offline, bakes
  result curves + rendered frames; HTML just displays with sliders snapping to the grid.
  + Most accurate, heaviest physics. - Sliders limited to precomputed values.
- **(B) Browser-live:** ship decimated geometry, do voxelize+PSF+metrics live in JS.
  + Fully continuous sliders. - JS perf risk; decimation fidelity risk.
- **(C) Hybrid (recommended, pending spike):** ship **decimated geometry for one
  representative subvolume** (axon centerlines, target < 5 MB) for the live interactive
  "microscope view"; **precompute population-scale yield/purity curves** across the
  NA x density grid for the quantitative panel. Best of both; honest about which numbers
  are live vs. precomputed.

### Data baking (192 MB cache -> browser asset)
- Reduce each axon's surface vertices to a **centerline polyline** (skeletonize or PCA
  /spline fit) or a heavy vertex decimation.
- Store as compact binary/JSON polylines; target whole-subvolume asset < 5 MB.
- Validate baked metrics against full-resolution Python ground truth (sanity gate).

### UI (vanilla HTML + 2D canvas, no heavy deps, single self-contained file)
- **Panel 1 - Microscope view:** live-rendered imaged plane of the labeled subset at
  current NA/density/depth; toggle ground-truth overlay; merged axons highlighted.
- **Sliders:** NA (drives both PSF axes), labeling density, imaging depth/plane;
  later optional: noise/dwell, edge-vs-center FOV.
- **Panel 2 - Yield vs purity:** the curve with the current operating point marked, and
  reference markers (Jeff 250/cell; proposal 10-200; cohort scenarios 1,100 / 4,400 /
  25,000).
- **Panel 3 - Readout:** labeled axons, resolvable fraction, estimated usable synapses,
  side-by-side with Jeff's sketch.

### Milestones
- **M0** spike (Part B) — decision gate.
- **M1** bake subvolume asset + static render in HTML.
- **M2** live sliders + confusion metric in-browser.
- **M3** population yield/purity curves + proposal annotations.
- **M4** polish; optional postsynaptic/dendrite extension (uses `data/dends/`).

### Risks / mitigations
- Decimation changes the metric -> validate baked vs full-res in Python (gate at M1).
- JS convolution too slow -> separable Gaussian on a downsampled grid; or fall to (A).
- Metric definition disputed -> spike surfaces it early; pick definition with you.

---

## D. Open decisions for you
1. **Confusion metric:** mixed-pixel mass vs. resolvable-axon-count — I recommend the
   combined "resolvable-axon yield." OK?
2. **Architecture:** A / B / C — I recommend **C (hybrid)**, confirmed after the spike.
3. **Scope now:** approve the **spike only** first, then decide on the full build from
   its results? (my recommendation) — or approve spike + full build together?
