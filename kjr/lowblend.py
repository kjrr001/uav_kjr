#!/usr/bin/env python3
"""Take the long wavelengths of one DSM and the short ones of another:
  out = detail + lowpass(shape - detail)          (Gaussian, sigma metres)
KJ-1: detail = own-BBA DSM (no image-locked rings, no east-strip steps), shape =
the ODM-pose DSM (v3), whose heights fit the 5 GCPs to 3.5 cm RMS - the own BBA
bends the block by +-10 cm against them (east strip, one-sided ties, GCP 1 in
3 images); with the 2-D grid correction the raw DSM bends +-30 cm. sigma 30 m: GCP
height RMS 5.2 cm (100 m: 7.5 cm) and the image-locked rings stay out (0.32 cm).
Usage: lowblend.py <detail dsm> <shape dsm> <out> [--sigma 30]
"""
import argparse
import numpy as np
from osgeo import gdal
from scipy import ndimage

gdal.UseExceptions()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("detail"); ap.add_argument("shape"); ap.add_argument("out")
    ap.add_argument("--sigma", type=float, default=30.0)
    a = ap.parse_args()
    d = gdal.Open(a.detail); gt = d.GetGeoTransform()
    zd = d.ReadAsArray().astype(np.float32)
    s = gdal.Warp("", a.shape, format="MEM", outputBounds=(gt[0], gt[3] + gt[5] * d.RasterYSize,
                  gt[0] + gt[1] * d.RasterXSize, gt[3]), width=d.RasterXSize, height=d.RasterYSize,
                  resampleAlg="bilinear", dstNodata=np.nan)
    zs = s.ReadAsArray().astype(np.float32)
    diff = zs - zd
    k = 16                                  # work on a 16x coarser grid (~2 m)
    H, W = diff.shape
    Hc, Wc = -(-H // k), -(-W // k)
    pad = np.full((Hc * k, Wc * k), np.nan, np.float32); pad[:H, :W] = diff
    blk = pad.reshape(Hc, k, Wc, k)
    m = np.isfinite(blk).sum((1, 3))
    c = np.where(m > 0, np.nansum(blk, (1, 3)) / np.maximum(m, 1), 0)
    w = (m > 0).astype(np.float32)
    sg = a.sigma / (gt[1] * k)
    low = ndimage.gaussian_filter(c * w, sg) / np.maximum(ndimage.gaussian_filter(w, sg), 1e-6)
    yy = (np.arange(H) + 0.5) / k - 0.5; xx = (np.arange(W) + 0.5) / k - 0.5
    corr = ndimage.map_coordinates(low, np.meshgrid(yy, xx, indexing="ij"), order=1, mode="nearest")
    out = np.where(np.isfinite(zd), zd + corr, np.nan).astype(np.float32)
    o = gdal.GetDriverByName("GTiff").Create(a.out, W, H, 1, gdal.GDT_Float32,
                                             ["TILED=YES", "COMPRESS=DEFLATE", "BIGTIFF=YES"])
    o.SetGeoTransform(gt); o.SetProjection(d.GetProjection()); o.GetRasterBand(1).WriteArray(out); o = None
    v = corr[np.isfinite(zd)]
    print(f"long-wave correction p1 {np.percentile(v, 1)*100:+.1f} p50 {np.median(v)*100:+.1f} p99 {np.percentile(v, 99)*100:+.1f} cm")


if __name__ == "__main__":
    main()
