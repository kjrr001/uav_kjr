#!/usr/bin/env python3
"""DSM height at each GCP minus the surveyed height (median of a 7x7 cell window).
GCP list: ODM gcp_list.txt (first line CRS, then x y z px py image [id]); one
value per GCP id (the column after the image name; the image name if absent).
Usage: gcp_check.py <gcp_list.txt> <dsm.tif> [<dsm.tif> ...]"""
import sys
import numpy as np
from osgeo import gdal

gdal.UseExceptions()
g = {}
for line in open(sys.argv[1]).read().splitlines()[1:]:
    p = line.split()
    if len(p) >= 6:
        g.setdefault(p[6] if len(p) > 6 else p[5], tuple(map(float, p[:3])))
for f in sys.argv[2:]:
    d = gdal.Open(f); gt = d.GetGeoTransform(); b = d.GetRasterBand(1)
    res = []
    for k, (x, y, h) in g.items():
        c = int((x - gt[0]) / gt[1]); r = int((y - gt[3]) / gt[5])
        if not (3 <= c < d.RasterXSize - 3 and 3 <= r < d.RasterYSize - 3):
            res.append(np.nan); continue
        w = b.ReadAsArray(c - 3, r - 3, 7, 7).astype(float); w[w < -1000] = np.nan
        res.append(np.nanmedian(w) - h)
    res = np.array(res)
    print(f, " ".join(f"{k}:{v * 100:+.1f}" for k, v in zip(g, res)),
          f"cm; RMS {np.sqrt(np.nanmean(res ** 2)) * 100:.1f} cm, mean {np.nanmean(res) * 100:+.1f} cm")
