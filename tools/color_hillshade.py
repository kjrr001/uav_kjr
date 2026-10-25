#!/usr/bin/env python3
"""Small colour hillshade of a DSM for a quick look (turbo colours x hillshade).
Heights are exaggerated (--exag) so cm-level ripples and steps show.
Usage: color_hillshade.py <dsm.tif> <out.jpg> [--res 0.76] [--exag 5] [--zmin z --zmax z]
Colour range default: 1st-99th percentile of the DSM."""
import argparse
import numpy as np
from osgeo import gdal
from matplotlib import cm
from PIL import Image

gdal.UseExceptions()
ap = argparse.ArgumentParser()
ap.add_argument("src"); ap.add_argument("dst")
ap.add_argument("--res", type=float, default=0.76, help="output pixel, m")
ap.add_argument("--exag", type=float, default=5.0)
ap.add_argument("--zmin", type=float); ap.add_argument("--zmax", type=float)
a = ap.parse_args()
ds = gdal.Warp("", a.src, format="MEM", xRes=a.res, yRes=a.res, resampleAlg="average")
z = ds.ReadAsArray().astype(np.float32); z[z < -1000] = np.nan
m = np.isfinite(z); zf = np.where(m, z, np.nanmedian(z))
lo = a.zmin if a.zmin is not None else np.nanpercentile(z, 1)
hi = a.zmax if a.zmax is not None else np.nanpercentile(z, 99)
gy, gx = np.gradient(zf * a.exag, a.res)
az, alt = np.radians(315), np.radians(45)
sl = np.arctan(np.hypot(gx, gy)); asp = np.arctan2(-gx, gy)
hs = np.clip(np.sin(alt) * np.cos(sl) + np.cos(alt) * np.sin(sl) * np.cos(az - asp), 0, 1)
c = cm.turbo(np.clip((zf - lo) / (hi - lo), 0, 1))[..., :3]
rgb = c * (0.35 + 0.65 * hs[..., None]); rgb[~m] = 1
Image.fromarray((rgb * 255).astype(np.uint8)).save(a.dst, quality=92)
print(f"{a.dst}: {z.shape[1]}x{z.shape[0]}, colours {lo:.2f}-{hi:.2f} m")
