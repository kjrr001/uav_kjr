#!/usr/bin/env python3
"""Re-make OpenSfM's undistorted images from the original JPEGs with a
better resampling kernel. OpenSfM brown model (normalised by max(w, h)):
  xn = x/z, r2 = xn^2 + yn^2, radial = 1 + k1 r2 + k2 r2^2 + k3 r2^3
  xd = xn radial + 2 p1 xn yn + p2 (r2 + 2 xn^2)
  yd = yn radial + p1 (r2 + 2 yn^2) + 2 p2 xn yn
  pixel = (f xd + c_x) * size + w/2 - 0.5     (same for y)
The undistorted camera is a pinhole with the same focal and a centred
principal point (opensfm/undistorted/reconstruction.json).
A camera from kjr/bba.py also carries kjr_dr: a radial pixel correction at the
brown position, u = u0 + dr(rho) (u0 - centre) / rho, and kjr_grid_x/y: a 2-D
bilinear pixel correction grid at the same position."""
import json
import os
import cv2
import numpy as np


def maps(cam, w, h):
    size = max(w, h)
    u, v = np.meshgrid(np.arange(w, dtype=np.float64), np.arange(h, dtype=np.float64))
    f = cam["focal_x"]
    xn = (u - w / 2 + 0.5) / size / f
    yn = (v - h / 2 + 0.5) / size / f
    r2 = xn * xn + yn * yn
    rad = 1 + cam["k1"] * r2 + cam["k2"] * r2 ** 2 + cam["k3"] * r2 ** 3
    xd = xn * rad + 2 * cam["p1"] * xn * yn + cam["p2"] * (r2 + 2 * xn * xn)
    yd = yn * rad + cam["p1"] * (r2 + 2 * yn * yn) + 2 * cam["p2"] * xn * yn
    mx = (f * xd + cam["c_x"]) * size + w / 2 - 0.5
    my = (f * yd + cam["c_y"]) * size + h / 2 - 0.5
    if "kjr_dr" in cam:          # radial residual correction from kjr/bba.py (px, piecewise linear in rho)
        dr = np.asarray(cam["kjr_dr"], float)
        du, dv = mx - (w / 2 - 0.5), my - (h / 2 - 0.5)
        rho = np.maximum(np.hypot(du, dv), 1e-6)
        x = rho / cam["kjr_dr_step"]
        i0 = np.minimum(x.astype(int), len(dr) - 2); a = x - i0
        d = ((1 - a) * dr[i0] + a * dr[i0 + 1]) / rho
        ux, uy = d * du, d * dv
        if "kjr_grid_x" in cam:      # 2-D bilinear grid correction (px), same brown position
            gy_, gx_ = cam["kjr_grid_shape"]; st = cam["kjr_grid_step"]
            for comp, key in ((0, "kjr_grid_x"), (1, "kjr_grid_y")):
                g = np.asarray(cam[key], np.float32).reshape(gy_, gx_)
                c = cv2.remap(g, (mx / st).astype(np.float32), (my / st).astype(np.float32),
                              cv2.INTER_LINEAR, borderMode=cv2.BORDER_REPLICATE)
                if comp == 0:
                    ux = ux + c
                else:
                    uy = uy + c
        mx, my = mx + ux, my + uy
    return mx.astype(np.float32), my.astype(np.float32)


def load_cam(run_dir):
    rec = json.load(open(os.path.join(run_dir, "opensfm/reconstruction.json")))[0]
    return next(iter(rec["cameras"].values()))


def main():
    """undistort.py <odm run> <original image dir> <out dir> [names...]"""
    import sys
    run, src_dir, out = sys.argv[1:4]
    names = sys.argv[4:] or sorted(os.listdir(os.path.join(run, "opensfm/undistorted/images")))
    os.makedirs(out, exist_ok=True)
    cam = load_cam(run); mx = my = None
    for n in names:
        dst = os.path.join(out, n)
        if os.path.exists(dst):
            continue
        img = cv2.imread(os.path.join(src_dir, n.replace(".tif", "")), cv2.IMREAD_UNCHANGED)
        if mx is None:
            mx, my = maps(cam, img.shape[1], img.shape[0])
        # Lanczos: OpenSfM's bilinear smooths the texture by an amount that
        # depends on the sub-pixel phase of the distortion shift, and that
        # phase pattern is concentric around the image centre -> matcher bias
        # in rings (KJ-1, 2026-10-25)
        cv2.imwrite(dst, cv2.remap(img, mx, my, cv2.INTER_LANCZOS4, borderMode=cv2.BORDER_CONSTANT))
        print(n, flush=True)


if __name__ == "__main__":
    main()
