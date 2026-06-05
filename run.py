#!/usr/bin/env python
"""One-command launcher for the synsim dashboard. Standard library only.

    python run.py                # static dashboard (no dependencies, works offline)
    python run.py --live         # full interactive dashboard (checks/installs deps)
    python run.py --live 8080    # ... on a custom port
    python run.py --no-browser   # don't auto-open a browser

STATIC mode serves the precomputed site in docs/ - it needs nothing but Python,
so a fresh clone runs immediately. LIVE mode runs the real 3D metric on any
geometry; it verifies the Python version, installs the scientific dependencies
into the current interpreter if they are missing, and needs the LICONN geometry
cache under data/ (not shipped in the repo).
"""
import argparse
import importlib
import subprocess
import sys
import webbrowser
from functools import partial
from http.server import ThreadingHTTPServer, SimpleHTTPRequestHandler
from pathlib import Path

ROOT = Path(__file__).resolve().parent
DOCS = ROOT / "docs"
LIVE_DEPS = ["numpy", "scipy", "skimage", "matplotlib", "PIL"]   # import names
PIP_DEPS = ["numpy", "scipy", "scikit-image", "matplotlib", "pillow"]
MIN_PY = (3, 9)        # dashboard + synsim run on 3.9; pyproject lists 3.10 aspirationally


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
        print(f"installing missing dependencies: {', '.join(PIP_DEPS)}")
        rc = subprocess.call([sys.executable, "-m", "pip", "install", *PIP_DEPS])
        if rc != 0:
            sys.exit("dependency install failed. Install manually:\n"
                     f"  {sys.executable} -m pip install {' '.join(PIP_DEPS)}")
        still = [m for m in LIVE_DEPS if importlib.util.find_spec(m) is None]
        if still:
            sys.exit(f"still missing after install: {', '.join(still)}")
    print("all live dependencies present.")


def serve_live(port, open_browser):
    ensure_live_deps()
    if not (ROOT / "data" / "axons").exists():
        print("\nWARNING: data/ LICONN geometry cache not found.")
        print("Live mode builds scenes from data/axons and data/dends, which are")
        print("not shipped in this repo. Without it, scene-building will fail.")
        print("Use the static dashboard instead:  python run.py\n")
    env = {**__import__("os").environ, "PYTHONPATH": str(ROOT)}
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
                    help="run the full interactive dashboard (installs deps; needs data/)")
    ap.add_argument("port", nargs="?", type=int, default=8000, help="port (default 8000)")
    ap.add_argument("--no-browser", action="store_true", help="do not open a browser")
    args = ap.parse_args()
    open_browser = not args.no_browser
    if args.live:
        serve_live(args.port, open_browser)
    else:
        serve_static(args.port, open_browser)


if __name__ == "__main__":
    main()
