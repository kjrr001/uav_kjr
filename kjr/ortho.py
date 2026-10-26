#!/usr/bin/env python3
"""ODM-KJR step 3: orthophoto from the fused DSM.

Every ortho pixel (grid of the ODM orthophoto, e.g. 6.5 cm) takes its height
from the fused DSM (bilinear), is projected into every undistorted image that
sees it, and is a weighted blend of those samples. Weight = (1 - r)^4 with r
the normalised distance from the image centre, so the most nadir image
dominates and image edges fade out without hard seams. Pixels outside the DSM
coverage are transparent, so ortho and DSM share one footprint.

Writes <out>/ortho.tif (RGBA, tiled, DEFLATE, EPSG:32638).
Usage: ortho.py <odm run> <dir with dsm_fused.tif> [--ortho ref_ortho.tif] [--img-dir dir]
"""
import argparse
import os
import sys

import cv2
import numpy as np
from osgeo import gdal

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from kjr.common import Run  # noqa: E402

gdal.UseExceptions()
TILE = 2048
POW = 4


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run"); ap.add_argument("out")
    ap.add_argument("--ortho", default=None)
    ap.add_argument("--img-dir", default=None, help="own undistorted images (kjr/undistort.py, Lanczos)")
    ap.add_argument("--threads", type=int, default=6)
    args = ap.parse_args()
    run = Run(args.run)
    run.img_dir = args.img_dir
    ref = gdal.Open(args.ortho or os.path.join(args.run, "odm_orthophoto/odm_orthophoto.tif"))
    gt = ref.GetGeoTransform(); W, H = ref.RasterXSize, ref.RasterYSize
    px = gt[1]

    dd = gdal.Open(os.path.join(args.out, "dsm_fused.tif"))
    dgt = dd.GetGeoTransform()
    dsm = dd.GetRasterBand(1).ReadAsArray()
    dmask = np.isfinite(dsm).astype(np.uint8)
    dsm_f = np.where(np.isfinite(dsm), dsm, np.nanmedian(dsm)).astype(np.float32)

    print(f"ortho grid {W}x{H} @ {px:.4f} m, DSM {dsm.shape[1]}x{dsm.shape[0]}; loading {len(run.names)} images", flush=True)
    imgs = {n: run.load_rgb(n) for n in run.names}
    # image footprint bboxes in UTM, from DSM median plane, for culling
    plane = run.ground_plane()
    bbox = {}
    for n in run.names:
        fp = run.footprint(n, plane) + run.offset
        bbox[n] = (fp[:, 0].min() - 20, fp[:, 0].max() + 20, fp[:, 1].min() - 20, fp[:, 1].max() + 20)

    drv = gdal.GetDriverByName("GTiff")
    o = drv.Create(os.path.join(args.out, "ortho.tif"), W, H, 4, gdal.GDT_Byte,
                   ["TILED=YES", "COMPRESS=DEFLATE", "BIGTIFF=YES", "PHOTOMETRIC=RGB", "ALPHA=YES"])
    o.SetGeoTransform(gt)
    o.SetProjection(ref.GetProjection())

    def tile(tr, tc):
        th, tw = min(TILE, H - tr), min(TILE, W - tc)
        cols, rows = np.meshgrid(np.arange(tc, tc + tw, dtype=np.float32) + 0.5,
                                 np.arange(tr, tr + th, dtype=np.float32) + 0.5)
        X = gt[0] + cols * px          # UTM
        Y = gt[3] - rows * px
        # DSM lookup (bilinear), DSM pixel centres at +0.5
        dc = ((X - dgt[0]) / dgt[1] - 0.5).astype(np.float32)
        dr = ((Y - dgt[3]) / dgt[5] - 0.5).astype(np.float32)
        Z = cv2.remap(dsm_f, dc, dr, cv2.INTER_LINEAR, borderMode=cv2.BORDER_REPLICATE)
        valid = cv2.remap(dmask, dc, dr, cv2.INTER_NEAREST, borderValue=0) > 0
        if not valid.any():
            return tr, tc, None, 0
        xa, xb, ya, yb = X.min(), X.max(), Y.min(), Y.max()
        acc = np.zeros((th, tw, 3), np.float32); wsum = np.zeros((th, tw), np.float32)
        Xl = (X - run.offset[0]).astype(np.float64); Yl = (Y - run.offset[1]).astype(np.float64)
        P = np.stack([Xl[valid], Yl[valid], Z[valid].astype(np.float64)], -1)
        used = 0
        for n in run.names:
            b = bbox[n]
            if b[1] < xa or b[0] > xb or b[3] < ya or b[2] > yb:
                continue
            cam = run.cams[n]
            uv, depth = cam.project(P)
            ru = np.abs(uv[:, 0] - cam.cx) / (cam.w / 2); rv = np.abs(uv[:, 1] - cam.cy) / (cam.h / 2)
            r = np.maximum(ru, rv)
            inside = (r < 0.995) & (depth > 0)
            if not inside.any():
                continue
            used += 1
            wv = np.where(inside, (1 - np.hypot(ru, rv) / np.sqrt(2)) ** POW, 0).astype(np.float32)
            mu = np.full((th, tw), -1, np.float32); mv = np.full((th, tw), -1, np.float32)
            mu[valid] = uv[:, 0]; mv[valid] = uv[:, 1]
            smp = cv2.remap(imgs[n], mu, mv, cv2.INTER_LINEAR, borderMode=cv2.BORDER_CONSTANT)
            wt = np.zeros((th, tw), np.float32); wt[valid] = wv
            acc += smp.astype(np.float32) * wt[..., None]; wsum += wt
        ok = valid & (wsum > 0)
        rgb = np.where(ok[..., None], acc / np.maximum(wsum, 1e-9)[..., None], 0)
        out = np.clip(rgb + 0.5, 0, 255).astype(np.uint8)
        return tr, tc, (out, (ok * 255).astype(np.uint8)), used

    # tiles in threads (cv2 / numpy release the GIL); writes stay in this thread
    from concurrent.futures import ThreadPoolExecutor
    jobs = [(tr, tc) for tr in range(0, H, TILE) for tc in range(0, W, TILE)]
    with ThreadPoolExecutor(args.threads) as ex:
        for tr, tc, res, used in ex.map(lambda j: tile(*j), jobs):
            if res is None:
                continue
            for k in range(3):
                o.GetRasterBand(k + 1).WriteArray(res[0][..., k], tc, tr)
            o.GetRasterBand(4).WriteArray(res[1], tc, tr)
            print(f"  tile {tr},{tc}: {used} images", flush=True)
    o = None


if __name__ == "__main__":
    main()
