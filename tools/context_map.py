#!/usr/bin/env python3
"""Context map: colour hillshade of a regional DEM (geographic, e.g. eoutput_hh.tif) with the
UAV block outlines (from their DSMs) as labelled red polygons; left = whole DEM, right = zoom.
Usage: context_map.py <regional dem.tif> <out.png> NAME=dsm.tif [NAME=dsm.tif ...]"""
import json, sys
import numpy as np
from osgeo import gdal
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource
from matplotlib.patches import Rectangle
gdal.UseExceptions()
sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
from footprint import footprint  # noqa: E402
dem, out = sys.argv[1:3]
ds = gdal.Open(dem); gt = ds.GetGeoTransform()
z = ds.GetRasterBand(1).ReadAsArray().astype(float)
nd = ds.GetRasterBand(1).GetNoDataValue()
if nd is not None: z[z == nd] = np.nan
z[z < -1000] = np.nan
H, W = z.shape
ext = [gt[0], gt[0] + W * gt[1], gt[3] + H * gt[5], gt[3]]
lat0 = (ext[2] + ext[3]) / 2
dx = gt[1] * 111320 * np.cos(np.radians(lat0)); dy = -gt[5] * 110950
polys = {kv.split("=", 1)[0]: footprint(kv.split("=", 1)[1], epsg=4326).tolist() for kv in sys.argv[3:]}
allc = np.vstack([np.array(c) for c in polys.values()])
m = 0.025
zx = [allc[:, 0].min() - m, allc[:, 0].max() + m]; zy = [allc[:, 1].min() - m, allc[:, 1].max() + m]

def shade(zz):
    lo, hi = np.nanpercentile(zz, [2, 98])
    ls = LightSource(315, 45)
    rgb = ls.shade(np.nan_to_num(zz, nan=lo), plt.cm.terrain, vert_exag=8, dx=dx, dy=dy,
                   vmin=lo, vmax=hi + (hi - lo) * 0.4, blend_mode="soft")
    rgb[~np.isfinite(zz)] = 1
    return rgb, lo, hi

fig, ax = plt.subplots(1, 2, figsize=(22, 10.5), gridspec_kw=dict(width_ratios=[1.25, 1]))
rgb, lo, hi = shade(z)
ax[0].imshow(rgb, extent=ext, aspect=1 / np.cos(np.radians(lat0)))
for b, c in polys.items():
    c = np.array(c); ax[0].fill(c[:, 0], c[:, 1], fc="none", ec="red", lw=1.2)
ax[0].add_patch(Rectangle((zx[0], zy[0]), zx[1] - zx[0], zy[1] - zy[0], fc="none", ec="red", lw=1.5, ls="--"))
ax[0].set_title("Eridu regional DEM (eoutput_hh) - UAV survey area", fontsize=14)
# zoom
c0 = int((zx[0] - gt[0]) / gt[1]); c1 = int((zx[1] - gt[0]) / gt[1]) + 1
r0 = int((zy[1] - gt[3]) / gt[5]); r1 = int((zy[0] - gt[3]) / gt[5]) + 1
sub = z[r0:r1, c0:c1]
rgb2, lo2, hi2 = shade(sub)
e2 = [gt[0] + c0 * gt[1], gt[0] + c1 * gt[1], gt[3] + r1 * gt[5], gt[3] + r0 * gt[5]]
ax[1].imshow(rgb2, extent=e2, aspect=1 / np.cos(np.radians(lat0)))
for b, c in polys.items():
    c = np.array(c)
    ax[1].fill(c[:, 0], c[:, 1], fc=(1, 0, 0, 0.12), ec="red", lw=2)
    cx, cy = c[:, 0].mean(), c[:, 1].max()
    ax[1].text(cx, cy + 0.002, b, color="red", fontsize=15, fontweight="bold", ha="center", va="bottom",
               bbox=dict(fc="white", ec="none", alpha=0.75, pad=1.5))
ax[1].set_xlim(e2[0], e2[1]); ax[1].set_ylim(e2[2], e2[3])
ax[1].set_title("UAV blocks (ODM-KJR DSM / ortho footprints)", fontsize=14)
for a in ax:
    a.set_xlabel("Longitude (°E)"); a.set_ylabel("Latitude (°N)"); a.grid(color="k", alpha=0.25, lw=0.5)
    a.ticklabel_format(useOffset=False)
for a, l, h in [(ax[0], lo, hi), (ax[1], lo2, hi2)]:
    sm = plt.cm.ScalarMappable(cmap=plt.cm.terrain, norm=plt.Normalize(l, h + (h - l) * 0.4))
    cb = fig.colorbar(sm, ax=a, shrink=0.6, pad=0.02); cb.set_label("Elevation (m)")
fig.savefig(out, dpi=110, bbox_inches="tight")
print("written", out)
