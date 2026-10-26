#!/usr/bin/env python3
"""Take the long wavelengths of one DSM and the short ones of another:
  out = detail + lowpass(shape - detail)          (Gaussian, sigma metres)
KJ-1: detail = own-BBA DSM (no image-locked rings, no east-strip steps), shape =
the ODM-pose DSM (v3), whose heights fit the 5 GCPs to 3.5 cm RMS - the own BBA
bends the block by +-10 cm against them (east strip, one-sided ties, GCP 1 in
3 images); with the 2-D grid correction the raw DSM bends +-30 cm. sigma 30 m: GCP
height RMS 5.2 cm (100 m: 7.5 cm) and the image-locked rings stay out (0.32 cm).
GCP datum (--gcp): the shape DSM carries ODM's few-cm vertical offset against the GCPs
(mean +5..+6 cm on KJ-2/3, JJ-2, -6 cm on JJ-1/3). One constant shift = -mean(DSM - GCP)
is applied when it is significant (|mean| > 2 standard errors). Leave-one-out over the
six blocks: an offset lowers the GCP RMS (KJ-3 7.1 -> 4.2, KJ-2 5.8 -> 3.3, JJ-2 7.4 ->
5.5 cm); a fitted plane over-fits 5 GCPs (JJ-3 9.9 -> 17.7 cm) and is not used.
Usage: lowblend.py <detail dsm> <shape dsm> <out> [--sigma 30] [--gcp gcp_list.txt]
"""
import argparse
import numpy as np
from osgeo import gdal
from scipy import ndimage

gdal.UseExceptions()


def gcp_residuals(path, z, gt):
    """DSM - GCP height per GCP id (median of 7x7 cells), as tools/gcp_check.py"""
    g = {}
    for line in open(path).read().splitlines()[1:]:
        p = line.split()
        if len(p) >= 6:
            g.setdefault(p[6] if len(p) > 6 else p[5], tuple(map(float, p[:3])))
    res = []
    for x, y, h in g.values():
        c = int((x - gt[0]) / gt[1]); r = int((y - gt[3]) / gt[5])
        if 3 <= c < z.shape[1] - 3 and 3 <= r < z.shape[0] - 3:
            v = np.nanmedian(z[r - 3:r + 4, c - 3:c + 4])
            if np.isfinite(v):
                res.append(v - h)
    return np.array(res)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("detail"); ap.add_argument("shape"); ap.add_argument("out")
    ap.add_argument("--sigma", type=float, default=30.0)
    ap.add_argument("--gcp", default=None, help="ODM gcp_list.txt: remove a significant mean offset")
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
    if a.gcp:
        r = gcp_residuals(a.gcp, out, gt)
        mu, se = r.mean(), r.std(ddof=1) / np.sqrt(len(r)) if len(r) > 1 else np.inf
        use = len(r) >= 3 and abs(mu) > 2 * se
        print(f"GCP: {len(r)} points, DSM - GCP mean {mu*100:+.1f} cm (se {se*100:.1f}) -> "
              + (f"shift {-mu*100:+.1f} cm" if use else "no shift (not significant)"))
        if use:
            out -= np.float32(mu)
    o = gdal.GetDriverByName("GTiff").Create(a.out, W, H, 1, gdal.GDT_Float32,
                                             ["TILED=YES", "COMPRESS=DEFLATE", "BIGTIFF=YES"])
    o.SetGeoTransform(gt); o.SetProjection(d.GetProjection()); o.GetRasterBand(1).WriteArray(out); o = None
    v = corr[np.isfinite(zd)]
    print(f"long-wave correction p1 {np.percentile(v, 1)*100:+.1f} p50 {np.median(v)*100:+.1f} p99 {np.percentile(v, 99)*100:+.1f} cm")


if __name__ == "__main__":
    main()
