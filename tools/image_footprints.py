#!/usr/bin/env python3
"""Per block: small ortho with the UAV block boundary (red) and the ground footprint of
every camera image (BBA poses, image border rays intersected with the DSM, lens distortion
ignored), labelled with the image number at the camera position.
Usage: image_footprints.py <name> <work dir> <products dir> <out.png> <small ortho.jpg>
  (products dir holds <name>_dsm.tif and <name>_ortho.tif; name = product prefix, e.g. KJ1_kjr)"""
import json, sys
import numpy as np
from osgeo import gdal
from scipy import ndimage
from scipy.spatial.transform import Rotation
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
gdal.UseExceptions()
sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
from footprint import footprint  # noqa: E402
B, run, prod, outpng, small_jpg = sys.argv[1:6]
rec = json.load(open(f"{run}/bba/odm/opensfm/reconstruction.json"))[0]
cam = next(iter(rec["cameras"].values())); w, h = cam["width"], cam["height"]; S = max(w, h)
f = cam["focal_x"] if "focal_x" in cam else cam["focal"]
off = np.array(list(map(float, open(f"{run}/bba/odm/odm_georeferencing/coords.txt").read().split("\n")[1].split()[:2])))
dsm = gdal.Open(f"{prod}/{B}_dsm.tif")
d2 = gdal.Translate("", dsm, format="MEM", xRes=1, yRes=1, resampleAlg="average")
gt = d2.GetGeoTransform(); Z = d2.ReadAsArray().astype(float); Z[~np.isfinite(Z) | (Z < -1000)] = np.nan
zmed = np.nanmedian(Z); Zf = np.where(np.isfinite(Z), Z, zmed)

def zat(x, y):
    c = (x - gt[0]) / gt[1]; r = (y - gt[3]) / gt[5]
    return ndimage.map_coordinates(Zf, [r, c], order=1, mode="nearest", cval=zmed)

n = 25
t = np.linspace(0, 1, n)
border = np.r_[np.c_[t * (w - 1), 0 * t], np.c_[(w - 1) + 0 * t, t * (h - 1)],
               np.c_[(1 - t) * (w - 1), (h - 1) + 0 * t], np.c_[0 * t, (1 - t) * (h - 1)]]
xn = np.c_[(border[:, 0] - (w - 1) / 2) / S / f, (border[:, 1] - (h - 1) / 2) / S / f, np.ones(len(border))]
fps = {}
for name, sh in sorted(rec["shots"].items()):
    R = Rotation.from_rotvec(sh["rotation"]).as_matrix(); tt = np.array(sh["translation"])
    C = -R.T @ tt; d = xn @ R            # R^T xn per row
    X = C + ((zmed - C[2]) / d[:, 2])[:, None] * d
    for _ in range(5):
        z = zat(X[:, 0] + off[0], X[:, 1] + off[1])
        X = C + ((z - C[2]) / d[:, 2])[:, None] * d
    fps[name] = (X[:, :2] + off, C[:2] + off)

# small ortho + its extent
o = gdal.Open(f"{prod}/{B}_ortho.tif"); og = o.GetGeoTransform()
import PIL.Image as Image
small = Image.open(small_jpg)
ext = [og[0], og[0] + og[1] * o.RasterXSize, og[3] + og[5] * o.RasterYSize, og[3]]
q = footprint(f"{prod}/{B}_dsm.tif")         # UTM, the DSM's CRS

allp = np.vstack([p for p, _ in fps.values()])
pad = 15
xl = [min(allp[:, 0].min(), ext[0]) - pad, max(allp[:, 0].max(), ext[1]) + pad]
yl = [min(allp[:, 1].min(), ext[2]) - pad, max(allp[:, 1].max(), ext[3]) + pad]
asp = (yl[1] - yl[0]) / (xl[1] - xl[0])
fig, ax = plt.subplots(figsize=(14, 14 * asp + 1) if asp < 1.6 else (10, 10 * asp + 1))
ax.imshow(np.asarray(small), extent=ext)
cmap = plt.cm.tab20
for k, (name, (p, c)) in enumerate(fps.items()):
    col = cmap(k % 20)
    ax.fill(p[:, 0], p[:, 1], fc="none", ec=col, lw=0.9, alpha=0.9)
    num = name.split("_")[-2] if name.count("_") >= 2 else name.split(".")[0].split("_")[-1]
    ax.plot(*c, "o", ms=2.5, color=col)
    ax.text(c[0], c[1], num, fontsize=6.5, color="k", ha="center", va="bottom",
            bbox=dict(fc="white", ec="none", alpha=0.6, pad=0.3))
ax.plot(q[:, 0], q[:, 1], color="red", lw=2.5)
ax.set_xlim(xl); ax.set_ylim(yl); ax.set_aspect("equal")
ax.ticklabel_format(useOffset=False, style="plain")
ax.set_xlabel("UTM 38N Easting (m)"); ax.set_ylabel("UTM 38N Northing (m)")
ax.set_title(f"{B}: ortho, UAV DSM/ortho boundary (red), {len(fps)} camera image footprints (BBA poses on DSM)", fontsize=11)
ax.grid(alpha=0.3)
fig.savefig(outpng, dpi=130, bbox_inches="tight"); print("written", outpng, len(fps))
