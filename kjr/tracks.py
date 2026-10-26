#!/usr/bin/env python3
"""ODM-KJR BBA step 1b: multi-view tie tracks from the pairwise ties of kjr/ties.py.

Why: the pairwise ties are 2-view (ALSC moves the second image's point onto the first
image's SIFT position, so a chain a-b-c over two pairs is two slightly different ground
points). kjr/bba.py grouped them only on the first image's SIFT position (star tracks),
so only 8-27 % of its points had 3+ views and the outer strips, which share few direct
pairs with the rest, were held weakly (JJ-3 west-edge stripes).

1. Every guided pair stores the SIFT keypoint index of each tie (ka, kb). The graph of
   (image, keypoint) nodes linked by ties is split into connected components = tracks.
   A component holding two keypoints of one image (a wrong link somewhere) is dropped.
2. Each track gets a reference: the member nearest its image centre, at its exact SIFT
   position. Every other member is re-measured from that reference by Gruen ALSC (the
   same kernel and gates as ties.py), batched per (reference image, image) pair, from
   its own SIFT position and the local affine of that image pair. So all observations
   of a track are of one ground point.
3. Members ALSC rejects are dropped; tracks with < 2 observations are dropped. Ties of
   pairs without keypoint indices (FLANN fallback) are kept as 2-view tracks.
Writes <out>/tracks.npz: names, cam, uv (x, y px), pt  (one row per observation).
Usage: tracks.py <stereo meta json> <original image dir> <ties out dir (has ties/, sift/)>
"""
import argparse, glob, json, os, sys, time
from collections import OrderedDict

import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import connected_components

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from kjr.ties import grey, sift_job          # noqa: E402

GATE = 2.0          # px, ALSC result vs its start (= ties.py)


class ImgCache:
    """standardised grey images (ALSC input), a few at a time"""
    def __init__(self, src, n=10):
        self.src, self.n, self.d = src, n, OrderedDict()

    def get(self, name):
        from tiepoints.tiepoints import _std
        if name not in self.d:
            if len(self.d) >= self.n:
                self.d.popitem(last=False)
            self.d[name] = np.nan_to_num(_std(grey(os.path.join(self.src, name))))
        self.d.move_to_end(name)
        return self.d[name]


def homography_affine(pl, pr):
    """per-point 2x2 Jacobian dpr/dpl of one homography fitted to the group (ALSC
    affine_init); tiepoints._local_affine's per-point least squares took minutes"""
    import cv2
    aff = np.tile([1.0, 0.0, 0.0, 1.0], (len(pl), 1))
    if len(pl) < 8:
        return aff
    H, _ = cv2.findHomography(pl, pr, cv2.LMEDS)
    if H is None:
        return aff
    x, y = pl[:, 0], pl[:, 1]
    w = H[2, 0] * x + H[2, 1] * y + H[2, 2]
    u = (H[0, 0] * x + H[0, 1] * y + H[0, 2]) / w
    v = (H[1, 0] * x + H[1, 1] * y + H[1, 2]) / w
    a11 = (H[0, 0] - u * H[2, 0]) / w; a12 = (H[0, 1] - u * H[2, 1]) / w
    a21 = (H[1, 0] - v * H[2, 0]) / w; a22 = (H[1, 1] - v * H[2, 1]) / w
    out = np.stack([a11, a12, a21, a22], 1)
    good = np.all(np.isfinite(out), 1) & (np.abs(out[:, 0] * out[:, 3] - out[:, 1] * out[:, 2] - 1) < 0.5)
    aff[good] = out[good]
    return aff


def alsc(L, R, pl, pr):
    """ALSC of pl (in L) into R from pr, as ties.refine; returns (q, ok)"""
    from tiepoints.tiepoints import alsc_kernel
    aff = homography_affine(pl, pr)
    for k in range(30):
        try:
            r = alsc_kernel("cuda").refine_batch(
                L, R, pl[:, 0], pl[:, 1], pr[:, 0], pr[:, 1],
                affine_init=aff, valid1=None, valid2=None,
                half_win=12, max_iter=60, conv_geo=0.01, conv_rad=0.01,
                model="translation", min_valid_fraction=0.9, affine_prior=0.0)
            break
        except RuntimeError:          # status -100: GPU memory held by another process
            time.sleep(5 + 10 * np.random.rand())
    else:
        raise RuntimeError("ALSC kept failing")
    q = np.stack([r.i2, r.j2], 1)
    ok = (np.asarray(r.status) == 0) & np.all(np.isfinite(q), 1) & (np.hypot(*(q - pr).T) <= GATE)
    return q, ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("meta"); ap.add_argument("src"); ap.add_argument("out")
    ap.add_argument("--sift-workers", type=int, default=4)
    a = ap.parse_args()
    t0 = time.time()
    tdir, cdir = os.path.join(a.out, "ties"), os.path.join(a.out, "sift")
    files = sorted(glob.glob(os.path.join(tdir, "*.npz")))
    pairs = [os.path.basename(f)[:-4].split("__") for f in files]
    names = sorted({n for p in pairs for n in p})
    idx = {n: k for k, n in enumerate(names)}
    os.makedirs(cdir, exist_ok=True)
    need = [n for n in names if not os.path.exists(os.path.join(cdir, n + ".d.npy"))]
    if need:                                   # the SIFT cache is removed after ties.py
        import multiprocessing as mp
        with mp.get_context("spawn").Pool(a.sift_workers) as pool:
            list(pool.imap_unordered(sift_job, [(n, a.src, cdir) for n in need]))
    xy = [np.load(os.path.join(cdir, n + ".xy.npy")) for n in names]
    off = np.concatenate([[0], np.cumsum([len(x) for x in xy])])

    # 1. keypoint graph
    ea, eb, loose = [], [], []
    for f, (x, y) in zip(files, pairs):
        with np.load(f) as d:
            if len(d["pa"]) < 30:
                continue
            if "ka" in d:
                ea.append(off[idx[x]] + d["ka"]); eb.append(off[idx[y]] + d["kb"])
            else:
                loose.append((idx[x], idx[y], d["pa"], d["pb"]))
    ea, eb = np.concatenate(ea), np.concatenate(eb)
    nodes = np.unique(np.r_[ea, eb])
    m = len(nodes)
    ia, ib = np.searchsorted(nodes, ea), np.searchsorted(nodes, eb)
    G = sparse.coo_matrix((np.ones(len(ia)), (ia, ib)), shape=(m, m))
    ncomp, lab = connected_components(G, directed=False)
    img = np.searchsorted(off, nodes, side="right") - 1
    kp = nodes - off[img]
    # components with two keypoints in one image
    key = lab.astype(np.int64) * len(names) + img
    u, c = np.unique(key, return_counts=True)
    bad = np.zeros(ncomp, bool); bad[u[c > 1] // len(names)] = True
    size = np.bincount(lab, minlength=ncomp)
    good = ~bad & (size >= 2)
    keep = good[lab]
    lab, img, kp = lab[keep], img[keep], kp[keep]
    hist = np.bincount(np.bincount(lab)[np.unique(lab)])
    print(f"keypoint graph: {m} nodes, {len(ea)} links, {ncomp} components, {bad.sum()} inconsistent; "
          f"kept {good.sum()} tracks, length histogram {dict((k, int(v)) for k, v in enumerate(hist) if v)}",
          flush=True)

    # 2. reference per track: the member nearest its image centre
    pos = np.array([xy[i][k] for i, k in zip(img, kp)], np.float64)
    import cv2
    h, w = cv2.imread(os.path.join(a.src, names[0]), cv2.IMREAD_REDUCED_GRAYSCALE_8).shape
    w, h = w * 8, h * 8
    rc = np.hypot(pos[:, 0] - w / 2, pos[:, 1] - h / 2)
    order = np.lexsort((rc, lab))
    lab, img, kp, pos = lab[order], img[order], kp[order], pos[order]
    first = np.r_[True, lab[1:] != lab[:-1]]
    tid = np.cumsum(first) - 1                          # 0..ntracks-1
    ref_row = np.nonzero(first)[0][tid]                 # row of each member's reference
    is_ref = first
    rimg, rpos = img[ref_row], pos[ref_row]

    # 3. ALSC from the reference into every other member, per (ref image, image)
    q = pos.copy(); ok = is_ref.copy()
    cache = ImgCache(a.src)
    oth = np.nonzero(~is_ref)[0]
    grp = rimg[oth] * len(names) + img[oth]
    o2 = np.argsort(grp, kind="stable"); oth, grp = oth[o2], grp[o2]
    cuts = np.r_[0, np.nonzero(np.diff(grp))[0] + 1, len(grp)]
    nmeas = 0
    tl = time.time()
    for gi, (s, e) in enumerate(zip(cuts[:-1], cuts[1:])):
        if gi % 50 == 0:
            print(f"  group {gi}/{len(cuts) - 1}, {time.time() - tl:.0f} s", flush=True)
        rows = oth[s:e]
        r_i, j_i = rimg[rows[0]], img[rows[0]]
        L, R = cache.get(names[r_i]), cache.get(names[j_i])
        qq, kk = alsc(L, R, rpos[rows], pos[rows])
        q[rows], ok[rows] = qq, kk
        nmeas += len(rows)
    print(f"ALSC from the references: {nmeas} members in {len(cuts) - 1} image pairs, "
          f"{ok[~is_ref].mean() * 100:.1f} % kept", flush=True)
    # tracks with >= 2 observations left
    cnt = np.bincount(tid, weights=ok, minlength=tid.max() + 1)
    sel = ok & (cnt[tid] >= 2)
    _, pt = np.unique(tid[sel], return_inverse=True)
    cam, uv = img[sel], q[sel]
    npt = pt.max() + 1
    L = np.bincount(np.bincount(pt))
    print(f"tracks: {npt}, observations {len(pt)}, views histogram "
          f"{dict((k, int(v)) for k, v in enumerate(L) if v)}", flush=True)
    # 4. fallback pairs as 2-view tracks
    for x, y, pa, pb in loose:
        n = len(pa)
        cam = np.r_[cam, np.full(n, x), np.full(n, y)]
        uv = np.r_[uv, pa, pb]
        pt = np.r_[pt, npt + np.arange(n), npt + np.arange(n)]
        npt += n
    np.savez(os.path.join(a.out, "tracks.npz"), names=np.array(names), cam=cam.astype(np.int32),
             uv=uv, pt=pt.astype(np.int64))
    print(f"wrote {a.out}/tracks.npz: {npt} points, {len(pt)} observations "
          f"({len(loose)} fallback pairs as 2-view), {time.time() - t0:.0f} s", flush=True)


if __name__ == "__main__":
    main()
