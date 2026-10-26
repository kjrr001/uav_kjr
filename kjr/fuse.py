#!/usr/bin/env python3
"""ODM-KJR step 2: fuse the per-pair height rasters into one DSM.

Per ground cell (grid of stereo_meta.json, 2 x ortho pixel):
  median over pairs -> keep values within max(3 * NMAD, 5 cm) of it ->
  weighted mean, weight = feather ramp from the pair raster's edge (pair
  edges are noisier, and a hard edge shows as a seam).
Pair bias: each pair's height raster is corrected by a plane (offset + tilt)
fitted to its difference from the previous fusion; the mean correction over
all pairs is forced to zero so the GCP datum is kept. Seams in the first
fusion came from these few-cm pair offsets / tilts.
Then: spikes against a 7x7 median (> 0.25 m) removed, holes filled from the
edges (GDAL FillNodata) inside the stereo coverage mask, nothing outside it.
Optional light smoothing (NaN-aware Gaussian, --sigma cells).

Writes <out>/dsm_fused.tif (smoothed, cell grid, NaN nodata), dsm_raw.tif (unsmoothed), npairs.tif,
       <out>/pair_planes.json (the corrections used).
Usage: fuse.py <stereo out dir> [--iter 2] [--sigma 0]
"""
import argparse
import glob
import json
import os

import cv2
import numpy as np
from osgeo import gdal
from scipy import ndimage

gdal.UseExceptions()
TILE = 2048
RAMP = 60          # cells (~8 m) from a pair raster's edge to full weight
ERODE = 30         # cells (~4 m) dropped at every pair raster edge: NMAD there is 10 -> 4 cm
                   # over the first 30 cells and single pairs swing -26..+37 cm (KJ-1)
B_REF = 127.0      # m; pairs with a shorter baseline get weight (B / B_REF)^2


BW = {}           # pair file -> baseline weight


def fill_holes(m):
    """= ndimage.binary_fill_holes (cross connectivity), by labelling the background once
    instead of iterated dilation (minutes on a 140 M cell grid)."""
    lab, _ = ndimage.label(~m)
    edge = np.unique(np.r_[lab[0], lab[-1], lab[:, 0], lab[:, -1]])
    return m | ~np.isin(lab, edge[edge > 0]) & (lab > 0)


def _diamond(r):
    y, x = np.mgrid[-r:r + 1, -r:r + 1]
    return (np.abs(x) + np.abs(y) <= r).astype(np.uint8)


def closing(m, it):
    """= ndimage.binary_closing(m, iterations=it): it x cross = one diamond of radius it
    (cv2, border = background as in scipy)."""
    k = _diamond(it); u = m.astype(np.uint8)
    d = cv2.dilate(u, k, borderType=cv2.BORDER_CONSTANT, borderValue=0)
    return cv2.erode(d, k, borderType=cv2.BORDER_CONSTANT, borderValue=0) > 0


def opening(m, it):
    k = _diamond(it); u = m.astype(np.uint8)
    e = cv2.erode(u, k, borderType=cv2.BORDER_CONSTANT, borderValue=0)
    return cv2.dilate(e, k, borderType=cv2.BORDER_CONSTANT, borderValue=0) > 0


def write_tif(path, arr, g, dtype=gdal.GDT_Float32):
    H, W = arr.shape
    ds = gdal.GetDriverByName("GTiff").Create(path, W, H, 1, dtype,
                                              ["TILED=YES", "COMPRESS=DEFLATE", "BIGTIFF=YES"])
    ds.SetGeoTransform((g["utm_x0"], g["cell"], 0, g["utm_y0"], 0, -g["cell"]))
    ds.SetProjection(g["wkt"])      # copied from the ODM ortho
    ds.GetRasterBand(1).WriteArray(arr)
    ds = None


def load_pair(f, plane, bw=1.0):
    with np.load(f) as d:
        z = d["z"].astype(np.float32); r0, c0 = int(d["r0"]), int(d["c0"])
        rad = d["rad"].astype(np.float32) if "rad" in d else None
    if plane is not None:
        a, b, c = plane
        h, w = z.shape
        z -= (a + b * (np.arange(w, dtype=np.float32) + c0)[None, :] / 1000
              + c * (np.arange(h, dtype=np.float32) + r0)[:, None] / 1000).astype(np.float32)
    fin = np.isfinite(z)
    # distance from the raster's outer edge only: internal match holes must not erode the data
    outer = fill_holes(closing(fin, 5))
    dist = ndimage.distance_transform_edt(outer)
    z[dist <= ERODE] = np.nan
    wgt = bw * np.clip((dist - ERODE) / RAMP, 0, 1).astype(np.float32) ** 2
    if rad is not None:
        # image corners: heights 5-10 cm low beyond radius ~1 (half width); keep them
        # for coverage but let any centre-of-image pair win where there is one
        wgt *= np.clip((1.15 - np.nan_to_num(rad, nan=2)) / 0.25, 0.05, 1)
    return r0, c0, z, wgt


def nanmedian_sorted(stack):
    """median along axis 0 ignoring NaN (sort puts NaN last)."""
    s = np.sort(stack, axis=0)
    n = np.isfinite(s).sum(0)
    lo = np.clip((n - 1) // 2, 0, None); hi = np.clip(n // 2, 0, None)
    a = np.take_along_axis(s, lo[None], 0)[0]; b = np.take_along_axis(s, hi[None], 0)[0]
    return np.where(n > 0, (a + b) / 2, np.nan), n


def fuse_once(win, planes, H, W):
    dsm = np.full((H, W), np.nan, np.float32)
    npr = np.zeros((H, W), np.uint8)
    cache = {}
    for tr in range(0, H, TILE):
        for tc in range(0, W, TILE):
            th, tw = min(TILE, H - tr), min(TILE, W - tc)
            hits = [x for x in win if x[1] < tr + th and x[1] + x[3] > tr and x[2] < tc + tw and x[2] + x[4] > tc]
            if not hits:
                continue
            stack = np.full((len(hits), th, tw), np.nan, np.float32)
            wst = np.zeros((len(hits), th, tw), np.float32)
            for k, (f, r0, c0, h, w) in enumerate(hits):
                if f not in cache:
                    if len(cache) > 60:
                        cache.pop(next(iter(cache)))
                    cache[f] = load_pair(f, planes.get(os.path.basename(f)), BW.get(os.path.basename(f), 1.0))
                _, _, z, wg = cache[f]
                ra, rb = max(tr, r0), min(tr + th, r0 + h)
                ca, cb = max(tc, c0), min(tc + tw, c0 + w)
                stack[k, ra - tr:rb - tr, ca - tc:cb - tc] = z[ra - r0:rb - r0, ca - c0:cb - c0]
                wst[k, ra - tr:rb - tr, ca - tc:cb - tc] = wg[ra - r0:rb - r0, ca - c0:cb - c0]
            med, n = nanmedian_sorted(stack)
            with np.errstate(all="ignore"):
                mad, _ = nanmedian_sorted(np.abs(stack - med))
                tol = np.maximum(3 * 1.4826 * np.nan_to_num(mad, nan=0), 0.05)
                keep = (np.abs(stack - med) <= tol)
                wk = np.where(keep, np.maximum(wst, 1e-3), 0)
                s = (np.where(keep, stack, 0) * wk).sum(0); c = wk.sum(0)
                dsm[tr:tr + th, tc:tc + tw] = np.where(c > 0, s / np.maximum(c, 1e-9), np.nan)
            npr[tr:tr + th, tc:tc + tw] = np.minimum(n, 255)
    return dsm, npr


def fit_planes(win, planes, dsm):
    """Robust plane (a + b x_km + c y_km, x/y in grid cells) of pair - dsm."""
    new = {}
    for f, r0, c0, h, w in win:
        key = os.path.basename(f)
        _, _, z, _ = load_pair(f, planes.get(key))
        dz = z - dsm[r0:r0 + h, c0:c0 + w]
        m = np.isfinite(dz)
        if m.sum() < 5000:
            new[key] = planes.get(key, [0, 0, 0]); continue
        yy, xx = np.nonzero(m)
        sel = slice(None, None, max(1, m.sum() // 200000))
        x, y, d = (xx[sel] + c0) / 1000, (yy[sel] + r0) / 1000, dz[m][sel].astype(float)
        keep = np.ones(len(d), bool)
        for _ in range(3):
            A = np.c_[np.ones(keep.sum()), x[keep], y[keep]]
            coef, *_ = np.linalg.lstsq(A, d[keep], rcond=None)
            r = d - (coef[0] + coef[1] * x + coef[2] * y)
            s = 1.4826 * np.median(np.abs(r[keep]))
            keep = np.abs(r) < 3 * s
        old = np.array(planes.get(key, [0, 0, 0]), float)
        new[key] = list(old + coef)
    # keep the datum: zero-mean correction over all pairs
    M = np.array(list(new.values())); M -= M.mean(0)
    return {k: list(v) for k, v in zip(new, M)}


def _prep(job):
    f, cdir = job
    key = os.path.basename(f)
    _, _, z, w = load_pair(f, None, BW.get(key, 1.0))
    np.save(os.path.join(cdir, key[:-4] + ".z.npy"), z)
    np.save(os.path.join(cdir, key[:-4] + ".w.npy"), w)
    return key


def _sorted_median(t, n):
    """t (K, ...) sorted with +inf for missing, n = count of finite: median as nanmedian_sorted."""
    import torch
    lo = ((n - 1) // 2).clamp(min=0); hi = (n // 2).clamp(min=0)
    a = torch.gather(t, 0, lo[None]).squeeze(0); b = torch.gather(t, 0, hi[None]).squeeze(0)
    return torch.where(n > 0, (a + b) / 2, torch.full_like(a, float("nan")))


def fuse_once_gpu(win, cdir, H, W, tile=1024):
    """fuse_once with the per-pair rasters prepared once (cache dir, memory-mapped) and the
    per-cell median / MAD / weighted mean on the GPU (torch). Same arithmetic as fuse_once."""
    import torch
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    dsm = np.full((H, W), np.nan, np.float32)
    npr = np.zeros((H, W), np.uint8)
    mm = {}
    for f, r0, c0, h, w in win:
        k = os.path.basename(f)[:-4]
        mm[f] = (np.load(os.path.join(cdir, k + ".z.npy"), mmap_mode="r"),
                 np.load(os.path.join(cdir, k + ".w.npy"), mmap_mode="r"))
    inf = float("inf")
    for tr in range(0, H, tile):
        for tc in range(0, W, tile):
            th, tw = min(tile, H - tr), min(tile, W - tc)
            hits = [x for x in win if x[1] < tr + th and x[1] + x[3] > tr and x[2] < tc + tw and x[2] + x[4] > tc]
            if not hits:
                continue
            stack = np.full((len(hits), th, tw), np.nan, np.float32)
            wst = np.zeros((len(hits), th, tw), np.float32)
            for k, (f, r0, c0, h, w) in enumerate(hits):
                z, wg = mm[f]
                ra, rb = max(tr, r0), min(tr + th, r0 + h)
                ca, cb = max(tc, c0), min(tc + tw, c0 + w)
                stack[k, ra - tr:rb - tr, ca - tc:cb - tc] = z[ra - r0:rb - r0, ca - c0:cb - c0]
                wst[k, ra - tr:rb - tr, ca - tc:cb - tc] = wg[ra - r0:rb - r0, ca - c0:cb - c0]
            S = torch.from_numpy(stack).to(dev); Wt = torch.from_numpy(wst).to(dev)
            fin = torch.isfinite(S)
            n = fin.sum(0)
            med = _sorted_median(torch.sort(torch.where(fin, S, inf), 0).values, n)
            dev_ = (S - med).abs()
            fd = torch.isfinite(dev_)
            mad = _sorted_median(torch.sort(torch.where(fd, dev_, inf), 0).values, fd.sum(0))
            tol = torch.clamp(3 * 1.4826 * torch.nan_to_num(mad, nan=0.0), min=0.05)
            keep = fd & (dev_ <= tol)
            wk = torch.where(keep, torch.clamp(Wt, min=1e-3), torch.zeros_like(Wt))
            sm = (torch.where(keep, S, torch.zeros_like(S)) * wk).sum(0); c = wk.sum(0)
            out = torch.where(c > 0, sm / torch.clamp(c, min=1e-9), torch.full_like(c, float("nan")))
            dsm[tr:tr + th, tc:tc + tw] = out.cpu().numpy()
            npr[tr:tr + th, tc:tc + tw] = n.clamp(max=255).to(torch.uint8).cpu().numpy()
    return dsm, npr


def median7(a, rows=256):
    """= ndimage.median_filter(a, size=7) (mode reflect = numpy symmetric), on the GPU in row bands."""
    import torch
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    H, W = a.shape
    p = np.pad(a, 3, mode="symmetric")
    out = np.empty_like(a)
    for r in range(0, H, rows):
        n = min(rows, H - r)
        t = torch.from_numpy(np.ascontiguousarray(p[r:r + n + 6])).to(dev)
        u = t.unfold(0, 7, 1).unfold(1, 7, 1).reshape(n, W, 49)
        out[r:r + n] = u.median(-1).values.cpu().numpy()
    return out


def finish(dsm, sigmas, npr):
    """spike removal + hole filling once; returns [(out, sigma) for each sigma], cov, spikes."""
    # coverage from where pairs contributed (geometry), not from surviving data:
    # noisy single-pair areas must not make a ragged edge
    cov = fill_holes(closing(npr > 0, 10))
    cov = opening(cov, 10)
    fin = np.isfinite(dsm)
    med7 = median7(np.where(fin, dsm, np.nanmedian(dsm)).astype(np.float32))
    spike = fin & (np.abs(dsm - med7) > 0.25)
    del med7
    dsm = dsm.copy(); dsm[spike] = np.nan
    H, W = dsm.shape
    mem = gdal.GetDriverByName("MEM").Create("", W, H, 1, gdal.GDT_Float32)
    b = mem.GetRasterBand(1); b.SetNoDataValue(-9999)
    b.WriteArray(np.where(np.isfinite(dsm), dsm, -9999).astype(np.float32))
    gdal.FillNodata(b, None, maxSearchDist=400, smoothingIterations=0)
    filled = b.ReadAsArray()
    base = np.where(cov & (filled > -9000), filled, np.nan).astype(np.float32)
    res = []
    for sigma in sigmas:
        out = base
        if sigma > 0:
            m = np.isfinite(base)
            num = ndimage.gaussian_filter(np.where(m, base, 0), sigma)
            den = ndimage.gaussian_filter(m.astype(np.float32), sigma)
            out = np.where(m, num / np.maximum(den, 1e-6), np.nan).astype(np.float32)
        res.append(out)
    return res, cov, spike.sum()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("--iter", type=int, default=0,
                    help="per-pair plane-correction iterations. 0 since 2026-10-25: on KJ-1 every pair"
                         " correction (this one, and a joint overlap adjustment of all pairs) made the"
                         " DSM worse against ODM ultra (east edge 6.6 -> 10.7 cm and more) - it created"
                         " the block steps it was meant to remove")
    ap.add_argument("--sigma", type=float, default=1.5,
                    help="final Gaussian, cells (1.5 = ~20 cm: looks like ODM ultra, keeps tracks/ditches)")
    ap.add_argument("--workers", type=int, default=12, help="processes preparing the pair rasters")
    ap.add_argument("--cpu", action="store_true", help="old single-process CPU fusion (reference)")
    args = ap.parse_args()
    out = args.out
    meta = json.load(open(os.path.join(out, "stereo_meta.json")))
    g = meta["grid"]
    for a, b, _ov, B in meta["pairs"]:
        BW[f"{a}__{b}.npz"] = min(1.0, (B / B_REF) ** 2)
    H, W = g["H"], g["W"]
    win = []
    for f in sorted(glob.glob(os.path.join(out, "pairs", "*.npz"))):
        with np.load(f) as d:
            win.append((f, int(d["r0"]), int(d["c0"]), *d["z"].shape))
    print(f"{len(win)} pair rasters, grid {W}x{H} @ {g['cell']:.4f} m", flush=True)

    planes = {}
    if args.iter == 0 and not args.cpu:
        import multiprocessing as mp
        import shutil
        cdir = os.path.join(out, ".fuse_cache")
        os.makedirs(cdir, exist_ok=True)
        with mp.get_context("fork").Pool(args.workers) as pool:
            for _ in pool.imap_unordered(_prep, [(x[0], cdir) for x in win], chunksize=4):
                pass
        dsm, npr = fuse_once_gpu(win, cdir, H, W)
        shutil.rmtree(cdir)
    for it in range(args.iter + 1) if (args.iter or args.cpu) else []:
        dsm, npr = fuse_once(win, planes, H, W)
        if it == args.iter:
            break
        planes = fit_planes(win, planes, dsm)
        P = np.array(list(planes.values()))
        print(f"  iter {it}: pair offset |p90| {np.percentile(np.abs(P[:, 0]), 90)*100:.1f} cm, "
              f"tilt p90 {np.percentile(np.hypot(P[:, 1], P[:, 2]) / g['cell'], 90):.2f} m/km", flush=True)
    json.dump(planes, open(os.path.join(out, "pair_planes.json"), "w"), indent=0)

    (raw, dsm), cov, nsp = finish(dsm, (0, args.sigma), npr)
    write_tif(os.path.join(out, "dsm_raw.tif"), raw, g)
    print(f"coverage {cov.sum() * g['cell']**2 / 1e6:.3f} km2, spikes removed {nsp}", flush=True)
    write_tif(os.path.join(out, "dsm_fused.tif"), dsm, g)
    write_tif(os.path.join(out, "npairs.tif"), npr, g, dtype=gdal.GDT_Byte)
    v = dsm[np.isfinite(dsm)]
    print(f"DSM z p1 {np.percentile(v, 1):.2f} p50 {np.median(v):.2f} p99 {np.percentile(v, 99):.2f}; "
          f"pairs/cell median {np.median(npr[cov]):.0f}", flush=True)


if __name__ == "__main__":
    main()
