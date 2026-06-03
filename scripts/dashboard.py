"""First-version dashboard for the sparse-axon imaging simulation.

A self-contained local web app (stdlib http.server, no extra deps). The browser
builds physics controls from `synsim.ui_spec()`; each change calls the backend,
which builds/caches the scene, runs the real 3D single-plane metric, and returns
the central optical section as raw data + bouton verdicts. The browser renders it
on a canvas, so display controls (LUT, min/max contrast) are instant and need no
server round-trip.

    python scripts/dashboard.py            # serve on http://localhost:8000
    python scripts/dashboard.py 8080       # custom port

Geometry knobs (FOV, bouton spacing/jitter) rebuild the scene (~15 s, cached);
NA / density / edges / purity / pixel re-score the cached scene (fast).
"""
import base64
import json
import sys
from http.server import ThreadingHTTPServer, BaseHTTPRequestHandler
from urllib.parse import urlparse, parse_qs

import numpy as np
import matplotlib

from synsim import SimParams, ui_spec, build_scene, evaluate, section_image

_SCENE_CACHE = {}


def _make_luts():
    luts = {}
    for name in ("gray", "viridis", "magma", "inferno", "hot", "cividis"):
        c = matplotlib.colormaps[name](np.linspace(0, 1, 256))[:, :3]
        luts[name] = (c * 255).astype(int).tolist()
    luts["green"] = [[0, i, 0] for i in range(256)]
    luts["red"] = [[i, 0, 0] for i in range(256)]
    return luts


LUTS = _make_luts()


def get_scene(params):
    key = (params.fov_xy_um, params.fov_z_um, params.bouton_spacing_um,
           params.bouton_jitter, int(round(params.region)))
    if key not in _SCENE_CACHE:
        _SCENE_CACHE[key] = build_scene(rng=np.random.default_rng(0), params=params)
    return _SCENE_CACHE[key]


def run(params):
    scene = get_scene(params)
    mask = np.random.default_rng(0).random(len(scene.axons)) < params.density
    r = evaluate(scene, mask, params)
    img, _ = section_image(scene, mask, params.na, params.pix_um)        # (nx, ny)
    mx = float(img.max())
    norm = (img / mx if mx > 0 else img).astype("<f4")                    # little-endian f32
    lo, fov = scene.origin, scene.fov

    pos, res, pur = r["positions"], r["resolved"], r["purity"]
    rel = (pos[:, :2] - lo[:2]) if len(pos) else np.zeros((0, 2))
    hist = (np.histogram(pur, bins=25, range=(0.0, 1.0))[0].tolist()
            if len(pur) else [0] * 25)

    return {
        "img": base64.b64encode(norm.tobytes(order="C")).decode(),
        "nx": int(img.shape[0]), "ny": int(img.shape[1]),
        "fov": [float(fov[0]), float(fov[1])],
        "edge_xy": params.edge_xy_um, "purity_thresh": params.purity_thresh,
        "bx": rel[:, 0].round(3).tolist(),
        "by": rel[:, 1].round(3).tolist(),
        "bres": [int(x) for x in res.tolist()],
        "hist": hist,
        "median_purity": (float(np.median(pur)) if len(pur) else None),
        "labeled_axons": r["labeled_axons"], "scored": r["scored"],
        "resolvable": r["resolvable"],
        "resolvable_frac": (None if r["scored"] == 0 else r["resolvable_frac"]),
        "synapses_per_um": scene.synapses_per_um, "n_axons": len(scene.axons),
    }


INDEX = r"""<!doctype html><html><head><meta charset=utf-8>
<title>synsim - sparse axon imaging</title>
<style>
 :root{--accent:#3b82f6;--ok:#16a34a;--bad:#dc2626;--bg:#0f1217;--panel:#171b22;--ink:#e8eaed;--mut:#8b94a3;--line:#2a313c}
 *{box-sizing:border-box}
 body{font:13px/1.4 system-ui,sans-serif;margin:0;display:flex;height:100vh;background:var(--bg);color:var(--ink)}
 #controls{width:330px;flex:none;padding:18px;overflow:auto;background:var(--panel);border-right:1px solid var(--line)}
 #controls::-webkit-scrollbar{width:8px} #controls::-webkit-scrollbar-thumb{background:var(--line);border-radius:4px}
 #view{flex:1;padding:22px;overflow:auto;display:flex;flex-direction:column;gap:16px;align-items:center}
 h1{font-size:15px;margin:0 0 4px;letter-spacing:.3px}
 .sub{color:var(--mut);font-size:11px;margin-bottom:14px}
 h2{font-size:10.5px;text-transform:uppercase;letter-spacing:.8px;color:var(--mut);margin:18px 0 8px;border-bottom:1px solid var(--line);padding-bottom:4px}
 .ctl{margin:9px 0} .ctl label{display:flex;justify-content:space-between;font-size:12.5px;margin-bottom:3px}
 .ctl input[type=range]{width:100%;accent-color:var(--accent);height:18px}
 .val{color:var(--accent);font-variant-numeric:tabular-nums;font-weight:600}
 .desc{color:var(--mut);font-size:10.5px;margin-top:1px}
 select{width:100%;background:#0f1217;color:var(--ink);border:1px solid var(--line);border-radius:6px;padding:5px}
 input[type=checkbox]{accent-color:var(--accent)}
 .card{background:var(--panel);border:1px solid var(--line);border-radius:12px;padding:14px}
 #stats{display:flex;gap:12px;width:100%;max-width:760px}
 .stat{flex:1;text-align:center;padding:10px 6px}
 .stat .n{font-size:24px;font-weight:700;font-variant-numeric:tabular-nums}
 .stat .l{font-size:10.5px;color:var(--mut);text-transform:uppercase;letter-spacing:.5px;margin-top:2px}
 #bigwrap .n{font-size:30px;color:var(--accent)}
 #cv{border-radius:10px;background:#000;image-rendering:pixelated;display:block}
 #hist{display:block;width:100%}
 .panel-title{font-size:11px;color:var(--mut);text-transform:uppercase;letter-spacing:.6px;margin-bottom:8px}
 #busy{position:fixed;top:12px;right:18px;background:var(--accent);color:#fff;padding:5px 12px;border-radius:20px;font-size:12px;display:none;box-shadow:0 2px 8px rgba(0,0,0,.4)}
</style></head><body>
<div id=controls>
 <h1>synsim</h1><div class=sub>sparse 2P axon imaging vs LICONN</div>
 <div id=display></div>
 <div id=form></div>
</div>
<div id=view>
 <div id=stats>
   <div class="card stat" id=bigwrap><div class=n id=big>-</div><div class=l>resolvable</div></div>
   <div class="card stat"><div class=n id=s_scored>-</div><div class=l>synapses imaged</div></div>
   <div class="card stat"><div class=n id=s_res>-</div><div class=l>resolvable</div></div>
   <div class="card stat"><div class=n id=s_lab>-</div><div class=l>labeled axons</div></div>
 </div>
 <div class=card><div class=panel-title>central optical section</div>
   <canvas id=cv width=720 height=720></canvas></div>
 <div class=card style="width:100%;max-width:760px"><div class=panel-title>per-bouton purity distribution (current settings)</div>
   <canvas id=hist width=720 height=190></canvas></div>
</div>
<div id=busy>computing...</div>
<script>
let spec=null, luts=null, vals={}, disp={lut:'gray',min:0,max:1,boutons:1}, last=null, timer=null;

async function init(){
 spec=await (await fetch('/api/spec')).json();
 luts=await (await fetch('/api/luts')).json();
 buildDisplay(); buildPhysics(); update();
}

function buildDisplay(){
 const d=document.getElementById('display');
 d.innerHTML='<h2>Display</h2>';
 let o=Object.keys(luts).map(n=>`<option ${n=='gray'?'selected':''}>${n}</option>`).join('');
 d.insertAdjacentHTML('beforeend',
  `<div class=ctl><label>LUT</label><select id=d_lut>${o}</select></div>
   <div class=ctl><label>min<span class=val id=v_min>0.00</span></label>
     <input type=range id=d_min min=0 max=1 step=0.01 value=0></div>
   <div class=ctl><label>max<span class=val id=v_max>1.00</span></label>
     <input type=range id=d_max min=0 max=1 step=0.01 value=1></div>
   <div class=ctl><label>show boutons<input type=checkbox id=d_bout checked></label></div>`);
 d_lut.onchange=e=>{disp.lut=e.target.value;redraw();};
 d_min.oninput=e=>{disp.min=+e.target.value;v_min.textContent=(+e.target.value).toFixed(2);redraw();};
 d_max.oninput=e=>{disp.max=+e.target.value;v_max.textContent=(+e.target.value).toFixed(2);redraw();};
 d_bout.onchange=e=>{disp.boutons=e.target.checked?1:0;redraw();};
}

function buildPhysics(){
 const groups={};
 for(const c of spec.controls){ vals[c.name]=c.value; (groups[c.group]??=[]).push(c); }
 const form=document.getElementById('form');
 for(const g in groups){ const h=document.createElement('h2'); h.textContent=g; form.appendChild(h);
  for(const c of groups[g]){
   const d=document.createElement('div'); d.className='ctl';
   d.innerHTML=`<label>${c.label}<span class=val id=p_${c.name}>${c.value}${c.unit?(' '+c.unit):''}</span></label>
    <input type=range min=${c.min} max=${c.max} step=${c.step} value=${c.value} id=i_${c.name}>
    <div class=desc>${c.desc}</div>`;
   form.appendChild(d);
   d.querySelector('input').addEventListener('input',e=>{
     vals[c.name]=parseFloat(e.target.value);
     document.getElementById('p_'+c.name).textContent=e.target.value+(c.unit?(' '+c.unit):'');
     clearTimeout(timer); timer=setTimeout(update,300);
   });
  }
 }
}

function b64f32(b64){ const s=atob(b64),n=s.length,u=new Uint8Array(n); for(let i=0;i<n;i++)u[i]=s.charCodeAt(i); return new Float32Array(u.buffer); }

async function update(){
 busy.style.display='block';
 const r=await (await fetch('/api/run?'+new URLSearchParams(vals))).json();
 if(r.error){ busy.textContent='error: '+r.error; return; }
 last={img:b64f32(r.img),nx:r.nx,ny:r.ny,fov:r.fov,edge:r.edge_xy,bx:r.bx,by:r.by,bres:r.bres,
       hist:r.hist,thr:r.purity_thresh,median:r.median_purity};
 big.textContent=r.resolvable_frac==null?'n/a':(100*r.resolvable_frac).toFixed(1)+'%';
 s_scored.textContent=r.scored; s_res.textContent=r.resolvable; s_lab.textContent=r.labeled_axons;
 redraw(); drawHist(); busy.style.display='none';
}

function redraw(){
 if(!last) return;
 const {img,nx,ny,fov,edge,bx,by,bres}=last, lut=luts[disp.lut];
 const cv=document.getElementById('cv'), c2=cv.getContext('2d');
 const off=document.createElement('canvas'); off.width=nx; off.height=ny;
 const oc=off.getContext('2d'), id=oc.createImageData(nx,ny);
 const span=(disp.max-disp.min)||1e-6;
 for(let iy=0;iy<ny;iy++)for(let ix=0;ix<nx;ix++){
   let t=(img[ix*ny+iy]-disp.min)/span; t=t<0?0:t>1?1:t;
   const col=lut[Math.min(255,Math.max(0,Math.floor(t*255)))];
   const p=(((ny-1-iy)*nx)+ix)*4;
   id.data[p]=col[0];id.data[p+1]=col[1];id.data[p+2]=col[2];id.data[p+3]=255;
 }
 oc.putImageData(id,0,0);
 c2.imageSmoothingEnabled=false; c2.clearRect(0,0,cv.width,cv.height);
 c2.drawImage(off,0,0,cv.width,cv.height);
 const ex=edge/fov[0]*cv.width, ey=edge/fov[1]*cv.height;
 c2.strokeStyle='rgba(11,191,255,.7)'; c2.setLineDash([5,4]); c2.lineWidth=1;
 c2.strokeRect(ex,ey,cv.width-2*ex,cv.height-2*ey); c2.setLineDash([]);
 if(disp.boutons) for(let k=0;k<bx.length;k++){
   const px=bx[k]/fov[0]*cv.width, py=cv.height-(by[k]/fov[1]*cv.height);
   c2.beginPath(); c2.arc(px,py,4.5,0,6.2832);
   c2.strokeStyle=bres[k]?'#39ff77':'#ff4d4d'; c2.lineWidth=1.5; c2.stroke();
 }
}

function drawHist(){
 if(!last||!last.hist) return;
 const cv=document.getElementById('hist'), ctx=cv.getContext('2d');
 const W=cv.width,H=cv.height,L=40,R=12,T=12,B=26;
 const pw=W-L-R, ph=H-T-B, n=last.hist.length, mx=Math.max(...last.hist,1);
 ctx.clearRect(0,0,W,H);
 ctx.fillStyle='#0f1217'; ctx.fillRect(L,T,pw,ph);          // plot area
 const bw=pw/n;
 for(let i=0;i<n;i++){
   const center=(i+0.5)/n, bh=(last.hist[i]/mx)*ph;
   ctx.fillStyle = center>=last.thr ? '#39c463' : '#e0564a';
   ctx.fillRect(L+i*bw, T+ph-bh, Math.max(1,bw-1), bh);
 }
 ctx.strokeStyle='#46505d'; ctx.lineWidth=1;                // axes
 ctx.beginPath(); ctx.moveTo(L,T); ctx.lineTo(L,T+ph); ctx.lineTo(L+pw,T+ph); ctx.stroke();
 const tx=L+last.thr*pw;                                    // threshold line
 ctx.strokeStyle='#f0f3f7'; ctx.setLineDash([4,3]);
 ctx.beginPath(); ctx.moveTo(tx,T); ctx.lineTo(tx,T+ph); ctx.stroke(); ctx.setLineDash([]);
 ctx.fillStyle='#c4ccd6'; ctx.font='11px system-ui';        // labels (light)
 ctx.textAlign='center';
 [0,0.25,0.5,0.75,1].forEach(v=>ctx.fillText(v.toFixed(2), L+v*pw, H-9));
 ctx.textAlign='right'; ctx.fillText(mx, L-4, T+10); ctx.fillText('0', L-4, T+ph);
 ctx.textAlign='left'; ctx.fillStyle='#f0f3f7';
 ctx.fillText('thr '+last.thr.toFixed(2)+(last.median!=null?'   median '+last.median.toFixed(2):''), tx+5, T+11);
 ctx.textAlign='left';
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
        elif u.path == "/api/luts":
            self._send(200, "application/json", json.dumps(LUTS).encode())
        elif u.path == "/api/run":
            q = {k: float(v[0]) for k, v in parse_qs(u.query).items()}
            try:
                self._send(200, "application/json", json.dumps(run(SimParams.from_dict(q))).encode())
            except Exception as e:
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
