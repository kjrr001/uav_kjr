"""Is the DSM ripple fixed to IMAGE coordinates (camera/lens) or to the GROUND?
High-pass the DSM (band 3-30 m), project every image's footprint into the image,
average the high-pass by image pixel bin over all images. Ground relief averages
out over 68 different footprints; a lens/camera error adds up.
Split images into two halves (even/odd) and correlate the two averages.
KJ-1 (2026-10-25): ODM / v3 0.84-0.92 corr (camera-model error), final v5 0.19.
Usage: diag_imgfixed.py <dsm.tif> <odm-like run (cameras)> <out.npy>   (+ out.png)"""
import sys, numpy as np
from osgeo import gdal
from scipy import ndimage
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from kjr.common import Run
dsm_p, odm = sys.argv[1], sys.argv[2]
ds = gdal.Open(dsm_p); gt = ds.GetGeoTransform()
f = 4  # 0.52 m cells
z = ds.ReadAsArray()[::1, ::1].astype(np.float32)
z = z[: z.shape[0] // f * f, : z.shape[1] // f * f].reshape(z.shape[0] // f, f, -1, f)
z[z < -1000] = np.nan
z = np.nanmean(z, (1, 3)); cell = gt[1] * f
m = np.isfinite(z); zf = np.where(m, z, np.nanmedian(z))
def g(s): 
    s /= cell; return ndimage.gaussian_filter(zf, s)
hp = g(1.5) - g(15.0); hp[~m] = np.nan; hp[np.abs(hp) > 0.15] = np.nan
run = Run(odm)
H, W = z.shape
yy, xx = np.mgrid[0:H, 0:W]
X = gt[0] + (xx + .5) * cell - run.offset[0]; Y = gt[3] - (yy + .5) * cell - run.offset[1]
P = np.c_[X.ravel(), Y.ravel(), zf.ravel()]; v = hp.ravel()
B = 64
c0 = run.cams[run.names[0]]
acc = [np.zeros((c0.h // B + 1, c0.w // B + 1)) for _ in range(2)]; cnt = [np.zeros_like(acc[0]) for _ in range(2)]
for k, n in enumerate(run.names):
    c = run.cams[n]
    uv, d = c.project(P)
    ok = (d > 0) & (uv[:, 0] >= 0) & (uv[:, 0] < c.w) & (uv[:, 1] >= 0) & (uv[:, 1] < c.h) & np.isfinite(v)
    i = (uv[ok, 1] // B).astype(int); j = (uv[ok, 0] // B).astype(int)
    np.add.at(acc[k % 2], (i, j), v[ok]); np.add.at(cnt[k % 2], (i, j), 1)
A = [a / np.maximum(c, 1) for a, c in zip(acc, cnt)]
ok = (cnt[0] > 50) & (cnt[1] > 50)
r = np.corrcoef(A[0][ok], A[1][ok])[0, 1]
avg = (A[0] + A[1]) / 2
print(f"hp DSM std {np.nanstd(hp)*100:.2f} cm; image-locked mean pattern std {avg[ok].std()*100:.2f} cm; even/odd corr {r:.2f}")
np.save(sys.argv[3], avg)
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 3, figsize=(18, 4.5))
for a, im, t in zip(ax, [A[0], A[1], avg], ["even images", "odd images", "all"]):
    h = a.imshow(np.where(ok, im, np.nan) * 100, cmap="RdBu_r", vmin=-3, vmax=3); a.set_title(t + " (cm, image coords)")
plt.colorbar(h, ax=ax); plt.savefig(sys.argv[3].replace(".npy", ".png"), dpi=90)
