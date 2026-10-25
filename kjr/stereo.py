#!/usr/bin/env python3
"""ODM-KJR step 1: pairwise dense stereo -> per-pair height rasters.

For every overlapping image pair of the ODM run (cameras unchanged):
  1. rectify both undistorted images to a common rotation (x = baseline);
     only the ground overlap of the pair is resampled;
  2. warp the right image by the disparity of the SfM ground plane (prior);
  3. census-SGM pyramid (UAV pfsgbm, soft census, MGM, CUDA) for the residual;
  4. triangulate, bin to the common ground grid (mean z per cell).
Writes <out>/pairs/<n1>__<n2>.npz = {r0, c0, z (float32, NaN), n (uint8)}.

Usage: stereo.py <odm run> <out dir> [--cell 0.13] [--only n1,n2] [--level 1] [--img-dir dir] [--meta-only]
"""
import argparse
import json
import os
import sys
import time

import cv2
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from kjr.common import Run, ray_plane, rect_rotation  # noqa: E402


def select_pairs(run, plane, min_overlap=0.2, min_base=30.0):
    fp = {n: run.footprint(n, plane).astype(np.float32) for n in run.names}
    area = {n: cv2.contourArea(fp[n]) for n in run.names}
    pairs = []
    for i, a in enumerate(run.names):
        for b in run.names[i + 1:]:
            B = np.linalg.norm(run.cams[a].C - run.cams[b].C)
            if B < min_base:
                continue
            inter, _ = cv2.intersectConvexConvex(cv2.convexHull(fp[a]), cv2.convexHull(fp[b]))
            ov = inter / min(area[a], area[b])
            if ov >= min_overlap:
                pairs.append((a, b, float(ov), float(B)))
    return pairs, fp


def grid_from_ortho(ortho_tif, offset, cell):
    from osgeo import gdal
    d = gdal.Open(ortho_tif)
    gt = d.GetGeoTransform()
    W = int(np.ceil(d.RasterXSize * gt[1] / cell))
    H = int(np.ceil(d.RasterYSize * -gt[5] / cell))
    # grid origin (top-left corner) in the run's local frame
    return dict(x0=gt[0] - offset[0], y0=gt[3] - offset[1], cell=cell, W=W, H=H,
                utm_x0=gt[0], utm_y0=gt[3], wkt=d.GetProjection())


EDGE_R = 9.0     # image radius / half width beyond which pixels are dropped. 1.0 made
                 # the block outline scalloped and lost 4 % coverage; corners are now
                 # down-weighted in the fusion instead (per-cell radius "rad" in the npz)


def cut_corners(img, cam):
    v, u = np.mgrid[0:img.shape[0], 0:img.shape[1]]
    img = img.copy()
    img[np.hypot(u - cam.cx, v - cam.cy) > EDGE_R * cam.w / 2] = np.nan
    return img


RECT_INTERP = cv2.INTER_LINEAR   # rectification resampling
PRIOR_ROUND = False              # round the plane-prior shift to whole pixels


def process_pair(run, a, b, plane, grid, matcher_kw, level=1, verbose=True):
    ca, cb = run.cams[a], run.cams[b]
    Rr = rect_rotation(ca, cb)
    B = np.linalg.norm(cb.C - ca.C)
    f = ca.f / level

    # ground overlap polygon -> rect canvases (same rows for both)
    fa = run.footprint(a, plane).astype(np.float32)
    fb = run.footprint(b, plane).astype(np.float32)
    _, poly = cv2.intersectConvexConvex(cv2.convexHull(fa), cv2.convexHull(fb))
    if poly is None or len(poly) < 3:
        return None
    poly = poly.reshape(-1, 2).astype(float)
    zpl = plane[0] + plane[1] * poly[:, 0] + plane[2] * poly[:, 1]
    P = np.vstack([np.c_[poly, zpl + dz] for dz in (-3.0, 3.0)])

    def rect_uv(C):
        Xr = (P - C) @ Rr.T
        return f * Xr[:, 0] / Xr[:, 2], f * Xr[:, 1] / Xr[:, 2]
    ua, va = rect_uv(ca.C)
    ub, _ = rect_uv(cb.C)
    pad = 8
    cxa, cxb, cyr = -ua.min() + pad, -ub.min() + pad, -va.min() + pad
    W = int(max(ua.max() - ua.min(), ub.max() - ub.min()) + 2 * pad)
    H = int(va.max() - va.min() + 2 * pad)
    Kr_a = np.array([[f, 0, cxa], [0, f, cyr], [0, 0, 1]])
    Kr_b = np.array([[f, 0, cxb], [0, f, cyr], [0, 0, 1]])

    def warp(cam, Kr, img):
        # rect pixel -> original pixel homography
        Hm = cam.K @ cam.R @ Rr.T @ np.linalg.inv(Kr)
        return cv2.warpPerspective(img, Hm, (W, H), flags=RECT_INTERP | cv2.WARP_INVERSE_MAP,
                                   borderMode=cv2.BORDER_CONSTANT, borderValue=np.nan)

    t0 = time.time()
    L = warp(ca, Kr_a, cut_corners(run.load_grey(a), ca))
    Rimg = warp(cb, Kr_b, cut_corners(run.load_grey(b), cb))

    # prior disparity from the plane: left canvas pixel -> plane -> right canvas column
    rows, cols = np.mgrid[0:H, 0:W].astype(np.float32)
    s = 4                                     # compute on a coarse lattice, smooth anyway
    rc, cc = rows[::s, ::s], cols[::s, ::s]
    d_rect = np.c_[(cc.ravel() - cxa) / f, (rc.ravel() - cyr) / f, np.ones(cc.size)]
    dw = d_rect @ Rr                           # world directions
    pa, pb, pc = plane
    t = (pa + pb * ca.C[0] + pc * ca.C[1] - ca.C[2]) / (dw[:, 2] - pb * dw[:, 0] - pc * dw[:, 1])
    X = ca.C + t[:, None] * dw
    Xb = (X - cb.C) @ Rr.T
    colb = f * Xb[:, 0] / Xb[:, 2] + cxb
    dprior = cv2.resize((colb - cc.ravel()).reshape(rc.shape).astype(np.float32), (W, H),
                        interpolation=cv2.INTER_LINEAR)
    if PRIOR_ROUND:
        # integer prior shift: no resampling of the right image here, so its
        # texture is not smoothed by a phase that varies across the canvas
        dprior = np.round(dprior)
    Rw = cv2.remap(Rimg, cols + dprior, rows, cv2.INTER_LINEAR, borderMode=cv2.BORDER_CONSTANT,
                   borderValue=np.nan)
    del Rimg

    # Single-level census-SGM. The plane prior already aligns the pair to
    # +-1 px, and the UAV pyramid driver drifts by +1..3 px on canvases larger
    # than ~1500 px (measured 2026-10-24 against phase correlation), so no pyramid.
    from matcher import fsgbm
    kw = dict(matcher_kw)
    srch, wsz = kw.pop("search_range"), kw.pop("window_size")
    Lf = np.where(np.isfinite(L), L, np.nanmean(L)).astype(np.float32)
    Rf = np.where(np.isfinite(Rw), Rw, np.nanmean(Rw)).astype(np.float32)
    ddx, _, score, ok = fsgbm.fsgbm_match(Lf, Rf, srch, wsz, **kw)
    del Lf, Rf
    # census window + cost box reach ~12 px: pixels that close to either image's
    # border matched against the mean fill, not image (seen as edge spikes)
    both = cv2.erode((np.isfinite(L) & np.isfinite(Rw)).astype(np.uint8), np.ones((25, 25), np.uint8)) > 0
    ok = ok & both & np.isfinite(ddx)
    t1 = time.time()

    # triangulate valid pixels
    colR = cols + dprior + ddx
    x1 = (cols - cxa) / f
    x2 = (colR - cxb) / f
    disp = x1 - x2
    ok &= disp > 1e-3
    Z = np.where(ok, B / np.where(ok, disp, 1), np.nan)
    Pr = np.stack([Z * x1, Z * (rows - cyr) / f, Z], -1)[ok]
    Xw = ca.C + Pr @ Rr                        # rect -> world
    del Pr

    # bin to the ground grid
    ci = np.floor((Xw[:, 0] - grid["x0"]) / grid["cell"]).astype(np.int64)
    ri = np.floor((grid["y0"] - Xw[:, 1]) / grid["cell"]).astype(np.int64)
    m = (ci >= 0) & (ci < grid["W"]) & (ri >= 0) & (ri < grid["H"])
    # image radius (/ half width) of each point, the larger of the two images:
    # heights fall 5-10 cm near the corners (KJ-1 radial profile), fuse.py down-weights them
    Xm = Xw[m]
    rad = np.zeros(len(Xm))
    for cam in (ca, cb):
        uv, _ = cam.project(Xm)
        rad = np.maximum(rad, np.hypot(uv[:, 0] - cam.cx, uv[:, 1] - cam.cy) / (cam.w / 2))
    ci, ri, z = ci[m], ri[m], Xw[m, 2]
    r0, c0 = ri.min(), ci.min()
    h, w = ri.max() - r0 + 1, ci.max() - c0 + 1
    lin = (ri - r0) * w + (ci - c0)
    n = np.bincount(lin, None, h * w)
    sz = np.bincount(lin, z, h * w)
    zz = np.where(n > 0, sz / np.maximum(n, 1), np.nan).reshape(h, w).astype(np.float32)
    rr = np.where(n > 0, np.bincount(lin, rad, h * w) / np.maximum(n, 1), np.nan).reshape(h, w).astype(np.float16)
    info = dict(canvas=(W, H), valid=float(ok.mean()), match_s=t1 - t0, B=B,
                z_med=float(np.nanmedian(zz)))
    if verbose:
        print(f"  {a} {b}: canvas {W}x{H}, valid {ok.mean()*100:.1f}%, B {B:.0f} m, "
              f"match {t1-t0:.1f}s, cells {np.isfinite(zz).sum()}", flush=True)
    return dict(r0=r0, c0=c0, z=zz, rad=rr, n=np.minimum(n, 255).reshape(h, w).astype(np.uint8)), info


# cost_box 1 + census 7: tuned 2026-10-25 by the agreement of two independent pairs
# (no shared image) over the same ground. cost_box 9 / census 9 had less fine noise
# but ~2 cm blob/ripple artifacts at 1.5-6 m that differ per pair (corr 0.41 ->
# 0.87, artifact 1.70 -> 0.46 cm); the fine noise is averaged by the fusion.
MATCHER = dict(search_range=8, window_size=9, cost="census", census_variant="soft",
               cost_box=1, census_win=7)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run"); ap.add_argument("out")
    ap.add_argument("--cell", type=float, default=0, help="default: 2 x ortho pixel")
    ap.add_argument("--ortho", default=None, help="ortho tif defining the grid extent")
    ap.add_argument("--only", default=None)
    ap.add_argument("--level", type=int, default=1)
    ap.add_argument("--max-pairs", type=int, default=0)
    ap.add_argument("--min-overlap", type=float, default=0.08)
    ap.add_argument("--img-dir", default=None, help="own undistorted images (kjr/undistort.py, Lanczos)")
    ap.add_argument("--meta-only", action="store_true", help="write stereo_meta.json (grid, pairs) and stop")
    args = ap.parse_args()

    run = Run(args.run)
    run.img_dir = args.img_dir
    plane = run.ground_plane()
    ortho = args.ortho or os.path.join(args.run, "odm_orthophoto/odm_orthophoto.tif")
    if not args.cell:
        from osgeo import gdal
        args.cell = 2 * gdal.Open(ortho).GetGeoTransform()[1]
    grid = grid_from_ortho(ortho, run.offset, args.cell)
    os.makedirs(os.path.join(args.out, "pairs"), exist_ok=True)
    pairs, _ = select_pairs(run, plane, min_overlap=args.min_overlap)
    if args.only:
        a, b = args.only.split(",")
        pairs = [p for p in pairs if {p[0], p[1]} == {a, b}] or [(a, b, 0, 0)]
    if args.max_pairs:
        pairs = pairs[:args.max_pairs]
    json.dump(dict(grid=grid, plane=list(plane), pairs=pairs, matcher=MATCHER, level=args.level),
              open(os.path.join(args.out, "stereo_meta.json"), "w"), indent=1, default=str)
    print(f"plane z = {plane[0]:.3f} + {plane[1]:.2e} x + {plane[2]:.2e} y; "
          f"grid {grid['W']}x{grid['H']} @ {args.cell} m; {len(pairs)} pairs", flush=True)
    if args.meta_only:
        return
    log = open(os.path.join(args.out, "pairs.log"), "a")
    for k, (a, b, ov, B) in enumerate(pairs):
        fn = os.path.join(args.out, "pairs", f"{a}__{b}.npz")
        if os.path.exists(fn):
            continue
        print(f"[{k+1}/{len(pairs)}] overlap {ov:.2f}", flush=True)
        res = process_pair(run, a, b, plane, grid, MATCHER, args.level)
        if res is None:
            continue
        arrs, info = res
        np.savez_compressed(fn, **arrs)
        log.write(json.dumps(dict(a=a, b=b, overlap=ov, **info)) + "\n"); log.flush()


if __name__ == "__main__":
    main()
