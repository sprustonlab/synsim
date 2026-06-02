"""Extract all LICONN axon meshes (expid82) to a compact, resumable cache.

Stores per-chunk .npz files with vertices in VOXEL coordinates (uint16) so the
whole set fits in ~0.4 GB. Faces are dropped (surface points suffice for PSF
voxelization; refetchable later if needed). Run with the isolated uv env:
    .venv-liconn/bin/python extract_axons.py
"""
import urllib.request
import time
from pathlib import Path

import numpy as np
from cloudvolume import CloudVolume

SRC = "precomputed://gs://liconn-public/ExPID82_1/segmentation/231030_agg_240123"
AXON_LIST = ("https://storage.googleapis.com/liconn-public/"
             "ExPID82_1/segmentation/231030_agg_240123/axons.txt")
OUT = Path("data/axons")
CHUNK = 200  # ids per fetch + per output file (resumable granularity)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    cv = CloudVolume(SRC, mip=0, use_https=True, progress=False)
    res_nm = np.asarray(cv.resolution, dtype=np.float64)  # [18,18,24]

    ids = [int(x) for x in urllib.request.urlopen(AXON_LIST, timeout=60).read().split()]
    print(f"{len(ids)} axon ids; res(nm)={res_nm.tolist()}; chunk={CHUNK}", flush=True)

    chunks = [ids[i:i + CHUNK] for i in range(0, len(ids), CHUNK)]
    t0 = time.time()
    n_done = 0
    n_verts = 0
    for ci, chunk in enumerate(chunks):
        out_f = OUT / f"chunk_{ci:04d}.npz"
        if out_f.exists():
            n_done += CHUNK
            continue
        meshes = cv.mesh.get(chunk)  # {id: Mesh}, missing ids simply absent
        verts_list, offsets, got_ids = [], [0], []
        for sid in chunk:
            m = meshes.get(sid)
            if m is None:
                continue
            v = np.asarray(m.vertices, dtype=np.float64) / res_nm  # nm -> voxel
            v = np.rint(v).astype(np.uint16)
            verts_list.append(v)
            got_ids.append(sid)
            offsets.append(offsets[-1] + v.shape[0])
        if verts_list:
            allv = np.concatenate(verts_list, axis=0)
        else:
            allv = np.zeros((0, 3), dtype=np.uint16)
        np.savez_compressed(
            out_f,
            verts=allv,
            offsets=np.asarray(offsets, dtype=np.int64),
            ids=np.asarray(got_ids, dtype=np.uint64),
        )
        n_done += len(chunk)
        n_verts += allv.shape[0]
        rate = n_done / max(time.time() - t0, 1e-6)
        eta = (len(ids) - n_done) / max(rate, 1e-6)
        print(f"[{ci+1}/{len(chunks)}] {n_done}/{len(ids)} axons "
              f"({rate:.0f}/s, ETA {eta/60:.1f} min), verts so far {n_verts:,}",
              flush=True)

    print(f"DONE in {(time.time()-t0)/60:.1f} min -> {OUT}/", flush=True)


if __name__ == "__main__":
    main()
