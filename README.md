# synsim — sparse 2P axon-imaging simulation vs LICONN

Simulates in-vivo two-photon imaging of **sparsely labeled axons** against real
LICONN ground-truth geometry, and asks how well individual presynaptic boutons
can be **resolved** as a function of numerical aperture, labeling density, and
the scoring criteria.

## Download and run locally (one command)

Clone the repo and run the launcher - standard library only, nothing to install
for the static dashboard:

```bash
git clone https://github.com/sprustonlab/synsim.git
cd synsim
python run.py                 # static dashboard at http://localhost:8000 (opens a browser)
```

`run.py` also drives the full interactive dashboard, checking the Python version
and installing the scientific dependencies if they are missing:

```bash
python run.py --live          # full dashboard; installs deps + downloads data
python run.py --live 8080     # ... on a custom port
python run.py --fetch-data    # only download the LICONN geometry cache, then exit
python run.py --no-browser    # don't auto-open a browser
```

The full dashboard needs the LICONN geometry cache under `data/`. If it is
missing, `--live` downloads it automatically from the **public**
`gs://liconn-public` bucket using the repo's `scripts/extract_*.py` (via
`cloud-volume`, installed on demand; resumable, ~11 min total). Pass `--no-fetch`
to skip this. The static mode needs none of it.

## Live demo (no install)

A fully static build of the dashboard is published via GitHub Pages:

**https://sprustonlab.github.io/synsim/**

It precomputes the optical-section image and per-bouton purity across a grid of
**NA x labeling density** (at the default field of view). The numerical aperture
and density sliders index that grid; the edge margins and the purity threshold
are applied live in the browser (a bouton's purity is independent of the
threshold and the margins, so no server is needed).

## Interactive dashboard (full, local)

The live dashboard runs the real 3D metric on any geometry, rebuilding scenes as
needed. It is a self-contained stdlib web app (no framework):

```bash
pip install -e .
python scripts/dashboard.py            # http://localhost:8000
python scripts/dashboard.py 8080       # custom port
```

Scene-building reads an extracted LICONN geometry cache under `data/`. It is not
committed to the repo, but `python run.py --fetch-data` (or `--live`) downloads
it from the public `gs://liconn-public` bucket via `scripts/extract_*.py`.

## Regenerating the static site

```bash
python scripts/export_static.py        # writes the grid + page into docs/
```

## Package layout

| Path | What |
|------|------|
| `synsim/` | imaging model, PSF, bouton placement, resolvability metric, params |
| `scripts/dashboard.py` | full interactive local dashboard (live compute) |
| `scripts/export_static.py` | bakes the static GitHub Pages build into `docs/` |
| `docs/` | the published static site (`index.html` + precomputed `cells/`) |
| `figures/` | analysis figures (NA / density / purity sweeps) |
