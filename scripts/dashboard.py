"""First-version dashboard for the sparse-axon imaging simulation.

A self-contained local web app (stdlib http.server, no extra deps). The browser
builds controls from `synsim.ui_spec()`; each change calls the backend, which
builds/caches the scene, runs the real 3D single-plane metric, and returns the
central optical section (boutons green=resolvable / red=confused) + the readout.

    python scripts/dashboard.py            # serve on http://localhost:8000
    python scripts/dashboard.py 8080       # custom port

Geometry knobs (FOV, bouton spacing/jitter) rebuild the scene (~15 s, cached);
NA / density / edges / purity / pixel re-score the cached scene (fast).
"""
import base64
import io
import json
import sys
from http.server import ThreadingHTTPServer, BaseHTTPRequestHandler
from urllib.parse import urlparse, parse_qs

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from synsim import SimParams, ui_spec, build_scene, evaluate, section_image

_SCENE_CACHE = {}


def get_scene(params):
    key = (params.fov_xy_um, params.fov_z_um, params.bouton_spacing_um, params.bouton_jitter)
    if key not in _SCENE_CACHE:
        _SCENE_CACHE[key] = build_scene(rng=np.random.default_rng(0), params=params)
    return _SCENE_CACHE[key]


def run(params):
    scene = get_scene(params)
    mask = np.random.default_rng(0).random(len(scene.axons)) < params.density
    r = evaluate(scene, mask, params)
    img, _ = section_image(scene, mask, params.na, params.pix_um)
    lo, fov = scene.origin, scene.fov

    fig, ax = plt.subplots(figsize=(6.4, 6.4))
    ax.imshow(img.T, origin="lower", cmap="gray", extent=[0, fov[0], 0, fov[1]])
    pos, res = r["positions"], r["resolved"]
    if len(pos):
        rel = pos - lo
        for good, c in [(res, "lime"), (~res, "red")]:
            ax.scatter(rel[good, 0], rel[good, 1], s=22, facecolor="none",
                       edgecolor=c, linewidths=0.9)
    # edge-exclusion box
    ex, ez = params.edge_xy_um, params.edge_z_um
    ax.add_patch(plt.Rectangle((ex, ex), fov[0] - 2 * ex, fov[1] - 2 * ex,
                               fill=False, edgecolor="deepskyblue", ls=":", lw=1))
    ax.set_xlim(0, fov[0]); ax.set_ylim(0, fov[1])
    ax.set_xlabel("x (um)"); ax.set_ylabel("y (um)")
    ax.set_title(f"central optical section (NA {params.na}, density {params.density:.2f})")
    buf = io.BytesIO()
    fig.tight_layout()
    fig.savefig(buf, format="png", dpi=100)
    plt.close(fig)
    b64 = base64.b64encode(buf.getvalue()).decode()

    return {
        "image": "data:image/png;base64," + b64,
        "labeled_axons": r["labeled_axons"],
        "scored": r["scored"],
        "resolvable": r["resolvable"],
        "resolvable_frac": (None if r["scored"] == 0 else r["resolvable_frac"]),
        "synapses_per_um": scene.synapses_per_um,
        "n_axons": len(scene.axons),
    }


INDEX = """<!doctype html><html><head><meta charset=utf-8>
<title>Sparse axon imaging - synsim</title>
<style>
 body{font:14px system-ui,sans-serif;margin:0;display:flex;height:100vh}
 #controls{width:340px;padding:16px;overflow:auto;border-right:1px solid #ddd;background:#fafafa}
 #view{flex:1;padding:16px;display:flex;flex-direction:column;gap:12px;align-items:center}
 h1{font-size:16px;margin:0 0 12px} h2{font-size:12px;text-transform:uppercase;color:#888;margin:14px 0 6px}
 .ctl{margin:8px 0} .ctl label{display:flex;justify-content:space-between;font-size:13px}
 .ctl input{width:100%} .val{color:#06c;font-variant-numeric:tabular-nums}
 .desc{color:#999;font-size:11px;margin-top:2px}
 #img{max-width:560px;width:100%;border:1px solid #ccc}
 #nums{display:grid;grid-template-columns:auto auto;gap:4px 18px;font-variant-numeric:tabular-nums}
 #nums .k{color:#666} #big{font-size:30px;font-weight:600;color:#06c}
 #busy{position:fixed;top:8px;right:12px;color:#c60;display:none}
</style></head><body>
<div id=controls><h1>synsim controls</h1><div id=form></div></div>
<div id=view>
 <div id=big>-</div>
 <img id=img>
 <div id=nums></div>
</div>
<div id=busy>computing...</div>
<script>
let spec=null, vals={}, timer=null;
async function init(){
 spec=await (await fetch('/api/spec')).json();
 const groups={};
 for(const c of spec.controls){ vals[c.name]=c.value; (groups[c.group]??=[]).push(c); }
 const form=document.getElementById('form');
 for(const g in groups){ const h=document.createElement('h2'); h.textContent=g; form.appendChild(h);
  for(const c of groups[g]){
   const d=document.createElement('div'); d.className='ctl';
   d.innerHTML=`<label>${c.label}<span class=val id=v_${c.name}>${c.value}${c.unit?(' '+c.unit):''}</span></label>
    <input type=range min=${c.min} max=${c.max} step=${c.step} value=${c.value} id=i_${c.name}>
    <div class=desc>${c.desc}</div>`;
   form.appendChild(d);
   d.querySelector('input').addEventListener('input',e=>{
     vals[c.name]=parseFloat(e.target.value);
     document.getElementById('v_'+c.name).textContent=e.target.value+(c.unit?(' '+c.unit):'');
     schedule();
   });
  }
 }
 update();
}
function schedule(){ clearTimeout(timer); timer=setTimeout(update,300); }
async function update(){
 document.getElementById('busy').style.display='block';
 const q=new URLSearchParams(vals).toString();
 const r=await (await fetch('/api/run?'+q)).json();
 document.getElementById('img').src=r.image;
 const pct=r.resolvable_frac==null?'n/a':(100*r.resolvable_frac).toFixed(1)+'%';
 document.getElementById('big').textContent='resolvable: '+pct;
 const f=x=>x==null?'-':x;
 document.getElementById('nums').innerHTML=
  `<div class=k>axons in FOV</div><div>${f(r.n_axons)}</div>
   <div class=k>labeled axons</div><div>${f(r.labeled_axons)}</div>
   <div class=k>boutons scored</div><div>${f(r.scored)}</div>
   <div class=k>resolvable</div><div>${f(r.resolvable)}</div>
   <div class=k>synapses / um</div><div>${r.synapses_per_um.toFixed(3)}</div>`;
 document.getElementById('busy').style.display='none';
}
init();
</script></body></html>"""


class Handler(BaseHTTPRequestHandler):
    def _send(self, code, ctype, body):
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self):
        u = urlparse(self.path)
        if u.path == "/":
            self._send(200, "text/html; charset=utf-8", INDEX.encode())
        elif u.path == "/api/spec":
            doc = {"defaults": SimParams().to_dict(), "controls": ui_spec()}
            self._send(200, "application/json", json.dumps(doc).encode())
        elif u.path == "/api/run":
            q = {k: float(v[0]) for k, v in parse_qs(u.query).items()}
            try:
                result = run(SimParams.from_dict(q))
                self._send(200, "application/json", json.dumps(result).encode())
            except Exception as e:  # surface errors to the browser
                self._send(500, "application/json", json.dumps({"error": str(e)}).encode())
        else:
            self._send(404, "text/plain", b"not found")

    def log_message(self, *a):
        pass


def main():
    port = int(sys.argv[1]) if len(sys.argv) > 1 else 8000
    print(f"synsim dashboard -> http://localhost:{port}  (first geometry change builds a scene, ~15 s)")
    ThreadingHTTPServer(("0.0.0.0", port), Handler).serve_forever()


if __name__ == "__main__":
    main()
