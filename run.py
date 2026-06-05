#!/usr/bin/env python
"""One-command launcher for the synsim dashboard. Standard library only.

    python run.py                # static dashboard (no dependencies, works offline)
    python run.py --live         # full interactive dashboard (installs deps + data)
    python run.py --live 8080    # ... on a custom port
    python run.py --fetch-data   # just download the LICONN geometry cache, then exit
    python run.py --no-browser   # don't auto-open a browser

STATIC mode serves the precomputed site in docs/ - it needs nothing but Python,
so a fresh clone runs immediately. LIVE mode runs the real 3D metric on any
geometry; it verifies the Python version, installs the scientific dependencies
if they are missing, and - if the LICONN geometry cache under data/ is absent -
downloads it from the public gs://liconn-public bucket via the repo's extract
scripts (resumable; ~11 min). Pass --no-fetch to skip that download.
"""
import argparse
import importlib
import os
import subprocess
import sys
import webbrowser
from functools import partial
from http.server import ThreadingHTTPServer, SimpleHTTPRequestHandler
from pathlib import Path

ROOT = Path(__file__).resolve().parent
DOCS = ROOT / "docs"
DATA = ROOT / "data"
LIVE_DEPS = ["numpy", "scipy", "skimage", "matplotlib", "PIL"]   # import names
PIP_DEPS = ["numpy", "scipy", "scikit-image", "matplotlib", "pillow"]
MIN_PY = (3, 9)        # dashboard + synsim run on 3.9; pyproject lists 3.10 aspirationally


def pip_install(packages):
    print(f"installing: {' '.join(packages)}")
    rc = subprocess.call([sys.executable, "-m", "pip", "install", *packages])
    if rc != 0:
        sys.exit("install failed. Run manually:\n"
                 f"  {sys.executable} -m pip install {' '.join(packages)}")


def data_present():
    """True if the LICONN geometry cache looks populated."""
    return (DATA / "axons").is_dir() and any((DATA / "axons").glob("chunk_*.npz"))


def fetch_data():
    """Download + cache the LICONN axon/dendrite geometry from the public bucket.

    Runs the repo's extract scripts, which pull meshes from the public
    gs://liconn-public bucket via cloud-volume. Resumable: existing chunks are
    skipped, so a partial download just continues.
    """
    print("LICONN geometry cache (data/) not found.")
    print("Downloading from the public bucket gs://liconn-public via cloud-volume")
    print("(~10 min for axons + ~1 min for dendrites; resumable).")
    if importlib.util.find_spec("cloudvolume") is None:
        pip_install(["cloud-volume"])
    env = {**os.environ, "PYTHONPATH": str(ROOT)}
    for script, label in (("extract_axons.py", "axons"), ("extract_dends.py", "dendrites")):
        print(f"\n--- fetching {label} ---", flush=True)
        rc = subprocess.call([sys.executable, str(ROOT / "scripts" / script)], env=env)
        if rc != 0:
            sys.exit(f"{script} failed (rc={rc}). You can re-run to resume.")
    print("\ngeometry cache ready in data/.")


def serve_static(port, open_browser):
    if not (DOCS / "index.html").exists():
        sys.exit(f"static site not found at {DOCS}/index.html (is this a full clone?)")
    handler = partial(SimpleHTTPRequestHandler, directory=str(DOCS))
    httpd = ThreadingHTTPServer(("0.0.0.0", port), handler)
    url = f"http://localhost:{port}/"
    print(f"synsim static dashboard -> {url}")
    print("(serving docs/ ; press Ctrl+C to stop)")
    if open_browser:
        webbrowser.open(url)
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        print("\nstopped.")


def ensure_live_deps():
    if sys.version_info < MIN_PY:
        sys.exit(f"live mode needs Python >= {MIN_PY[0]}.{MIN_PY[1]} "
                 f"(this is {sys.version.split()[0]}). Try: python run.py  (static mode)")
    missing = [m for m in LIVE_DEPS if importlib.util.find_spec(m) is None]
    if missing:
        pip_install(PIP_DEPS)
        still = [m for m in LIVE_DEPS if importlib.util.find_spec(m) is None]
        if still:
            sys.exit(f"still missing after install: {', '.join(still)}")
    print("all live dependencies present.")


def serve_live(port, open_browser, fetch=True):
    ensure_live_deps()
    if not data_present():
        if fetch:
            fetch_data()
        else:
            print("\nWARNING: data/ LICONN geometry cache not found and --no-fetch set.")
            print("Scene-building will fail. Run:  python run.py --fetch-data\n")
    env = {**os.environ, "PYTHONPATH": str(ROOT)}
    if open_browser:
        # dashboard prints its URL; open after a beat so the server is up.
        import threading
        import time

        def _open():
            time.sleep(2.0)
            webbrowser.open(f"http://localhost:{port}/")
        threading.Thread(target=_open, daemon=True).start()
    sys.exit(subprocess.call(
        [sys.executable, str(ROOT / "scripts" / "dashboard.py"), str(port)], env=env))


def main():
    ap = argparse.ArgumentParser(description="Launch the synsim dashboard.")
    ap.add_argument("--live", action="store_true",
                    help="run the full interactive dashboard (installs deps; fetches data/)")
    ap.add_argument("--fetch-data", action="store_true",
                    help="just download the LICONN geometry cache into data/, then exit")
    ap.add_argument("--no-fetch", action="store_true",
                    help="with --live, do not auto-download data/ if missing")
    ap.add_argument("port", nargs="?", type=int, default=8000, help="port (default 8000)")
    ap.add_argument("--no-browser", action="store_true", help="do not open a browser")
    args = ap.parse_args()
    if args.fetch_data:
        if data_present():
            print("data/ already present; nothing to do.")
        else:
            fetch_data()
        return
    open_browser = not args.no_browser
    if args.live:
        serve_live(args.port, open_browser, fetch=not args.no_fetch)
    else:
        serve_static(args.port, open_browser)


if __name__ == "__main__":
    main()
