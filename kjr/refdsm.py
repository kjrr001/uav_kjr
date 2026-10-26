#!/usr/bin/env python3
"""ODM-KJR without ODM: a reference DSM from the adjusted tie points (bba.py points.npy).

Replaces ODM's DSM where kjr used it: the centre of the +-8 px stereo search (stereo.py
--prior-dsm, smoothed to 1 m there anyway) and the long-wave shape of lowblend.py.
Median height per CELL (default 2 m) cell, cells with < 3 points empty, holes filled
(GDAL FillNodata), Gaussian sigma 1 cell. UTM, written into <run>/odm_dem/dsm.tif.
Usage: refdsm.py <bba out dir> <run dir (coords.txt, kjr_exif.json EPSG)> [--cell 2]
"""
import argparse, json, os

import numpy as np
from osgeo import gdal, osr
from scipy import ndimage

gdal.UseExceptions()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("bba"); ap.add_argument("run")
    ap.add_argument("--cell", type=float, default=2.0)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    X = np.load(os.path.join(a.bba, "points.npy")).astype(np.float64)
    with open(os.path.join(a.run, "odm_georeferencing/coords.txt")) as f:
        f.readline(); off = np.array(list(map(float, f.readline().split()[:2])))
    ej = os.path.join(a.run, "kjr_exif.json")
    epsg = json.load(open(ej))["epsg"] if os.path.exists(ej) else 32638
    # robust: drop points far from a local median (wrong ties)
    x, y, z = X[:, 0] + off[0], X[:, 1] + off[1], X[:, 2]
    c = a.cell
    x0, y1 = np.floor(x.min() / c) * c, np.ceil(y.max() / c) * c
    W, H = int(np.ceil((x.max() - x0) / c)) + 1, int(np.ceil((y1 - y.min()) / c)) + 1
    ci = ((x - x0) / c).astype(int); ri = ((y1 - y) / c).astype(int)
    lin = ri * W + ci
    o = np.lexsort((z, lin)); lin, zs = lin[o], z[o]
    u, st, cnt = np.unique(lin, return_index=True, return_counts=True)
    med = zs[st + cnt // 2]
    dem = np.full(H * W, np.nan, np.float32)
    ok = cnt >= 3
    dem[u[ok]] = med[ok]
    dem = dem.reshape(H, W)
    # spikes against a 5x5 median
    fin = np.isfinite(dem)
    m5 = ndimage.median_filter(np.where(fin, dem, np.nanmedian(dem)), 5)
    dem[fin & (np.abs(dem - m5) > 0.5)] = np.nan
    mem = gdal.GetDriverByName("MEM").Create("", W, H, 1, gdal.GDT_Float32)
    b = mem.GetRasterBand(1); b.SetNoDataValue(-9999)
    b.WriteArray(np.where(np.isfinite(dem), dem, -9999))
    gdal.FillNodata(b, None, maxSearchDist=50, smoothingIterations=0)
    f = b.ReadAsArray()
    m = f > -9000
    num = ndimage.gaussian_filter(np.where(m, f, 0), 1.0); den = ndimage.gaussian_filter(m.astype(np.float32), 1.0)
    out = np.where(m, num / np.maximum(den, 1e-6), np.nan).astype(np.float32)
    dst = a.out or os.path.join(a.run, "odm_dem/dsm.tif")
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    d = gdal.GetDriverByName("GTiff").Create(dst, W, H, 1, gdal.GDT_Float32, ["COMPRESS=DEFLATE", "TILED=YES"])
    d.SetGeoTransform((x0, c, 0, y1, 0, -c))
    sr = osr.SpatialReference(); sr.ImportFromEPSG(epsg); d.SetProjection(sr.ExportToWkt())
    d.GetRasterBand(1).SetNoDataValue(np.nan); d.GetRasterBand(1).WriteArray(out); d = None
    print(f"{len(X)} tie points -> {W}x{H} @ {c} m, {ok.sum()} cells with >= 3 points "
          f"({ok.sum() / max(m.sum(), 1) * 100:.0f} % of the filled area) -> {dst}")


if __name__ == "__main__":
    main()
