#!/usr/bin/env python3
"""ODM-KJR BBA step 1: sub-pixel tie points on the ORIGINAL (distorted) images.

SIFT -> guided GPU matching -> MAGSAC++ (F, 1 px) -> Gruen ALSC (cuda, UAV library).
  1. SIFT once per image (cached in <out>/sift, uint8 descriptors - OpenCV's are integers).
  2. Per pair, guided exact matching on the GPU: a homography from the 8000 (else 30000)
     strongest keypoints of each image, then each of 8 x 6 tiles of the first image is matched
     (exact 2-NN, ratio 0.8) only against the second image's keypoints inside the tile's
     predicted footprint +-120 px. Replaces FLANN over all ~500 k x 500 k keypoints of the
     UAV library's match_pair: KJ-2 pairs 21-29 s -> 0.4-0.9 s, and ~2x the MAGSAC inliers
     (16 k -> 30 k), as FLANN's approximate neighbours miss many true ones.
     If the homography fails the pair falls back to tiepoints.match_pair (FLANN).
  3. MAGSAC++ and ALSC exactly as tiepoints.match_pair.
Pairs: the overlapping pairs of a stereo run (stereo_meta.json), including cross-strip pairs.
Each pair keeps its own 2-view ties (ALSC moves only the second image's point; multi-view
tracks are grouped on the first image's unmoved SIFT position, kjr/bba.py).
Writes <out>/ties/<a>__<b>.npz with pa, pb (n,2) px (x = column, y = row).
Usage: ties.py <stereo meta json> <original image dir> <out> [--workers 2] [--flann]
"""
import argparse, json, os, shutil, sys, time
import multiprocessing as mp
import cv2
import numpy as np

SIFT = dict(nfeatures=0, contrastThreshold=0.012, edgeThreshold=8)   # = UAV default_args (sift)
RATIO = 0.8
NCOARSE, GRID, MARGIN = (8000, 30000), (8, 6), 120   # coarse tries: 8000 strongest, else 30000
                                                      # (KJ-3 cross-strip pairs with 10 % overlap)


def grey(path):
    im = cv2.imread(path, cv2.IMREAD_COLOR).astype(np.float32)
    return im @ np.array([0.114, 0.587, 0.299], np.float32)


def sift_job(t):
    name, src, cdir = t
    base = os.path.join(cdir, name)
    if os.path.exists(base + ".d.npy"):
        return name, None
    from features.feature_match import norm01
    t0 = time.time()
    A = grey(os.path.join(src, name))
    k, d = cv2.SIFT_create(**SIFT).detectAndCompute((norm01(A) * 255).astype(np.uint8), None)
    np.save(base + ".xy.npy", np.float32([p.pt for p in k]))
    np.save(base + ".rsp.npy", np.float32([p.response for p in k]))
    np.save(base + ".d.npy", d.astype(np.uint8))      # integer valued, <= 255
    return name, (len(k), time.time() - t0)


def _nn2(L, R, Lh, Rh, chunk=4096):
    """exact 2-NN of the rows of L among the rows of R: fp16 cosine top-3, then exact fp32."""
    import torch
    I, D = [], []
    for s in range(0, L.shape[0], chunk):
        top = (Lh[s:s + chunk] @ Rh.T).topk(min(3, R.shape[0]), dim=1).indices
        d = ((R[top] - L[s:s + chunk, None, :]) ** 2).sum(-1)
        d, o = d.sort(1)
        I.append(torch.gather(top, 1, o)); D.append(d)
    return torch.cat(I), torch.cat(D)


def guided_match(xyL, dL, rL, xyR, dR, rR, shape):
    import torch
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    L = torch.from_numpy(dL).to(dev).float(); R = torch.from_numpy(dR).to(dev).float()
    Lh = (L / L.norm(dim=1, keepdim=True).clamp(min=1e-6)).half()
    Rh = (R / R.norm(dim=1, keepdim=True).clamp(min=1e-6)).half()
    H = None
    for nc in NCOARSE:
        a = np.argsort(-rL)[:nc]; b = np.argsort(-rR)[:nc]
        ta, tb = torch.from_numpy(a).to(dev), torch.from_numpy(b).to(dev)
        idx, d = _nn2(L[ta], R[tb], Lh[ta], Rh[tb])
        ok = (d[:, 0].sqrt() < RATIO * d[:, 1].sqrt()).cpu().numpy()
        if ok.sum() < 30:
            continue
        j = idx[:, 0].cpu().numpy()
        H, m = cv2.findHomography(xyL[a[ok]], xyR[b[j[ok]]], cv2.RANSAC, 25.0)
        if H is not None and m.sum() >= 30:
            break
        H = None
    if H is None:
        return None
    h, w = shape
    tx = np.minimum((xyL[:, 0] / w * GRID[0]).astype(int), GRID[0] - 1)
    ty = np.minimum((xyL[:, 1] / h * GRID[1]).astype(int), GRID[1] - 1)
    qi, ti = [], []
    for gy in range(GRID[1]):
        for gx in range(GRID[0]):
            sel = np.nonzero((tx == gx) & (ty == gy))[0]
            if not len(sel):
                continue
            c = np.float32([[gx, gy], [gx + 1, gy], [gx, gy + 1], [gx + 1, gy + 1]]) * [w / GRID[0], h / GRID[1]]
            p = cv2.perspectiveTransform(c[None], H)[0]
            x0, y0 = p.min(0) - MARGIN; x1, y1 = p.max(0) + MARGIN
            cand = np.nonzero((xyR[:, 0] >= x0) & (xyR[:, 0] <= x1) & (xyR[:, 1] >= y0) & (xyR[:, 1] <= y1))[0]
            if len(cand) < 2:
                continue
            ts, tc = torch.from_numpy(sel).to(dev), torch.from_numpy(cand).to(dev)
            idx, d = _nn2(L[ts], R[tc], Lh[ts], Rh[tc])
            ok = d[:, 0].sqrt() < RATIO * d[:, 1].sqrt()
            qi.append(ts[ok].cpu().numpy()); ti.append(tc[idx[ok, 0]].cpu().numpy())
    if not qi:
        return None
    qi = np.concatenate(qi); ti = np.concatenate(ti)
    u, cnt = np.unique(ti, return_counts=True)        # one keypoint taken by several: ambiguous
    keep = ~np.isin(ti, u[cnt > 1])
    return xyL[qi[keep]].astype(np.float64), xyR[ti[keep]].astype(np.float64)


def refine(A, B, pl, pr, st):
    """MAGSAC++ and ALSC, as tiepoints.match_pair."""
    from tiepoints.tiepoints import _std, _local_affine, alsc_kernel
    st["candidates"] = len(pl)
    F, m = cv2.findFundamentalMat(pl.astype(np.float32), pr.astype(np.float32),
                                  cv2.USAC_MAGSAC, 1.0, 0.9999, 20000)
    if m is None:
        return np.zeros((0, 2)), np.zeros((0, 2))
    m = m.ravel().astype(bool)
    pl, pr = pl[m], pr[m]
    st["ransac_F"] = int(m.sum())
    L, R = _std(A), _std(B)
    aff = _local_affine(pl, pr)
    r = alsc_kernel("cuda").refine_batch(
        np.nan_to_num(L), np.nan_to_num(R), pl[:, 0], pl[:, 1], pr[:, 0], pr[:, 1],
        affine_init=aff, valid1=np.isfinite(L), valid2=np.isfinite(R),
        half_win=12, max_iter=60, conv_geo=0.01, conv_rad=0.01,
        model="translation", min_valid_fraction=0.9, affine_prior=0.0)
    q = np.stack([r.i2, r.j2], 1)
    good = (np.asarray(r.status) == 0) & np.all(np.isfinite(q), 1) & (np.hypot(*(q - pr).T) <= 2.0)
    st["final"] = int(good.sum())
    return pl[good], q[good]


def job(t):
    a, b, src, out, cdir, flann = t
    dst = os.path.join(out, f"{a}__{b}.npz")
    if os.path.exists(dst):
        return a, b, None
    t0 = time.time()
    A, B = grey(os.path.join(src, a)), grey(os.path.join(src, b))
    for k in range(20):
        try:
            st = {}
            res = None
            if not flann:
                ld = lambda n, e: np.load(os.path.join(cdir, n + e))
                res = guided_match(ld(a, ".xy.npy"), ld(a, ".d.npy"), ld(a, ".rsp.npy"),
                                   ld(b, ".xy.npy"), ld(b, ".d.npy"), ld(b, ".rsp.npy"), A.shape)
                import torch
                torch.cuda.empty_cache()      # give the GPU back to the ALSC kernel (own allocator)
            if res is None:
                from tiepoints.tiepoints import match_pair
                tp = match_pair(A, B, backend="sift", alsc="cuda", max_points=None)
                pa, pb, st = tp.pi, tp.pj, dict(tp.stats, fallback=1)
            else:
                pa, pb = refine(A, B, res[0], res[1], st)
            break
        except RuntimeError:          # ALSC status -100 / CUDA OOM: another worker holds the GPU
            import torch
            torch.cuda.empty_cache()
            time.sleep(5 + 10 * np.random.rand())
    else:
        raise RuntimeError(f"{a} {b}: kept failing")
    np.savez(dst, pa=np.asarray(pa, np.float64), pb=np.asarray(pb, np.float64))
    return a, b, (st.get("ransac_F", 0), st.get("final", 0), time.time() - t0, st.get("fallback", 0))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("meta"); ap.add_argument("src"); ap.add_argument("out")
    ap.add_argument("--workers", type=int, default=2, help="pair processes sharing the GPU")
    ap.add_argument("--sift-workers", type=int, default=4)
    ap.add_argument("--flann", action="store_true", help="old path: tiepoints.match_pair (FLANN) for every pair")
    ap.add_argument("--keep-sift", action="store_true", help="keep the SIFT cache (~80 MB per image)")
    a = ap.parse_args()
    tdir = os.path.join(a.out, "ties"); cdir = os.path.join(a.out, "sift")
    os.makedirs(tdir, exist_ok=True); os.makedirs(cdir, exist_ok=True)
    pairs = [(p[0].replace(".tif", ""), p[1].replace(".tif", "")) for p in json.load(open(a.meta))["pairs"]]
    todo = [p for p in pairs if not os.path.exists(os.path.join(tdir, f"{p[0]}__{p[1]}.npz"))]
    ctx = mp.get_context("spawn")                      # fork breaks the CUDA ALSC
    if todo and not a.flann:
        names = sorted({n for p in todo for n in p})
        with ctx.Pool(a.sift_workers) as pool:
            for n, s in pool.imap_unordered(sift_job, [(n, a.src, cdir) for n in names]):
                if s:
                    print(f"sift {n}: {s[0]} keypoints {s[1]:.1f}s", flush=True)
    jobs = [(x, y, a.src, tdir, cdir, a.flann) for x, y in pairs]
    with ctx.Pool(a.workers) as pool:
        for x, y, s in pool.imap_unordered(job, jobs):
            print(x, y, "skip" if s is None else
                  f"F {s[0]} kept {s[1]} {s[2]:.1f}s" + (" (FLANN fallback)" if s[3] else ""), flush=True)
    if not a.keep_sift:
        shutil.rmtree(cdir, ignore_errors=True)


if __name__ == "__main__":
    main()
