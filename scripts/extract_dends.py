"""Extract all LICONN dendrite meshes (expid82) to a compact cache.

Same scheme as extract_axons.py (vertices in voxel coords, uint16, faces
dropped). Run with the isolated uv env:
    .venv-liconn/bin/python extract_dends.py
"""
import urllib.request
import time
from pathlib import Path

import numpy as np
from cloudvolume import CloudVolume

SRC = "precomputed://gs://liconn-public/ExPID82_1/segmentation/231030_agg_240123"
DEND_LIST = ("https://storage.googleapis.com/liconn-public/"
             "ExPID82_1/segmentation/231030_agg_240123/dendrites.txt")
OUT = Path(__file__).resolve().parent.parent / "data" / "dends"
CHUNK = 200


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    cv = CloudVolume(SRC, mip=0, use_https=True, progress=False)
    res_nm = np.asarray(cv.resolution, dtype=np.float64)
    ids = [int(x) for x in urllib.request.urlopen(DEND_LIST, timeout=60).read().split()]
    print(f"{len(ids)} dendrite ids; res(nm)={res_nm.tolist()}", flush=True)

    chunks = [ids[i:i + CHUNK] for i in range(0, len(ids), CHUNK)]
    t0 = time.time()
    n_done = n_verts = 0
    for ci, chunk in enumerate(chunks):
        out_f = OUT / f"chunk_{ci:04d}.npz"
        if out_f.exists():
            n_done += CHUNK
            continue
        meshes = cv.mesh.get(chunk)
        verts_list, offsets, got_ids = [], [0], []
        for sid in chunk:
            m = meshes.get(sid)
            if m is None:
                continue
            v = np.rint(np.asarray(m.vertices, dtype=np.float64) / res_nm).astype(np.uint16)
            verts_list.append(v)
            got_ids.append(sid)
            offsets.append(offsets[-1] + v.shape[0])
        allv = np.concatenate(verts_list, axis=0) if verts_list else np.zeros((0, 3), np.uint16)
        np.savez_compressed(out_f, verts=allv,
                            offsets=np.asarray(offsets, np.int64),
                            ids=np.asarray(got_ids, np.uint64))
        n_done += len(chunk)
        n_verts += allv.shape[0]
        print(f"[{ci+1}/{len(chunks)}] {n_done}/{len(ids)} dends, verts {n_verts:,}", flush=True)
    print(f"DONE in {(time.time()-t0)/60:.1f} min -> {OUT}/", flush=True)


if __name__ == "__main__":
    main()
