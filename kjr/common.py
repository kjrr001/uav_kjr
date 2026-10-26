"""ODM-KJR common pieces: cameras from the ODM/OpenSfM undistorted
reconstruction, image loading, ground plane, rectifying rotation.

Frame: OpenSfM world of the ODM run = UTM minus the offset in
odm_georeferencing/coords.txt (x east, y north, z up, metres).
Camera: X_cam = R X + t, pixel = f * X_cam[:2] / X_cam[2] + (w/2 - 0.5, h/2 - 0.5)
(OpenSfM perspective model on the undistorted images, focal normalised by max(w, h)).
"""
import json
import os
from dataclasses import dataclass

import cv2
import numpy as np


@dataclass
class Cam:
    name: str
    R: np.ndarray      # world -> camera
    C: np.ndarray      # centre, world
    f: float           # px
    cx: float
    cy: float
    w: int
    h: int

    @property
    def K(self):
        return np.array([[self.f, 0, self.cx], [0, self.f, self.cy], [0, 0, 1.0]])

    def project(self, X):
        """X (n,3) world -> (n,2) px, depth (n,)"""
        Xc = (X - self.C) @ self.R.T
        z = Xc[:, 2]
        return np.c_[self.f * Xc[:, 0] / z + self.cx, self.f * Xc[:, 1] / z + self.cy], z


class Run:
    """An ODM project directory (read only)."""

    def __init__(self, run_dir):
        self.dir = run_dir
        rec = json.load(open(os.path.join(run_dir, "opensfm/undistorted/reconstruction.json")))[0]
        self.cams = {}
        for name, s in rec["shots"].items():
            cam = rec["cameras"][s["camera"]]
            w, h = cam["width"], cam["height"]
            R, _ = cv2.Rodrigues(np.array(s["rotation"], float))
            t = np.array(s["translation"], float)
            self.cams[name] = Cam(name, R, -R.T @ t, cam["focal"] * max(w, h),
                                  w / 2 - 0.5, h / 2 - 0.5, w, h)
        self.names = sorted(self.cams)
        self.sparse = np.array([p["coordinates"] for p in rec["points"].values()], float)
        with open(os.path.join(run_dir, "odm_georeferencing/coords.txt")) as f:
            zone = f.readline().split()[-1]            # "WGS84 UTM 38N"
            self.offset = np.array(list(map(float, f.readline().split()[:2])))
        self.epsg = (32600 if zone[-1].upper() == "N" else 32700) + int(zone[:-1])

    img_dir = None    # override: own undistorted images (kjr/undistort.py), same names

    def image_path(self, name):
        return os.path.join(self.img_dir or os.path.join(self.dir, "opensfm/undistorted/images"), name)

    def load_rgb(self, name):
        img = cv2.imread(self.image_path(name), cv2.IMREAD_UNCHANGED)
        return cv2.cvtColor(img, cv2.COLOR_BGR2RGB)

    def load_grey(self, name):
        rgb = self.load_rgb(name).astype(np.float32)
        return rgb @ np.array([0.299, 0.587, 0.114], np.float32)

    def ground_plane(self):
        """Robust plane z = a + b x + c y through the sparse SfM points."""
        P = self.sparse
        m = np.ones(len(P), bool)
        for _ in range(5):
            A = np.c_[np.ones(m.sum()), P[m, 0], P[m, 1]]
            coef, *_ = np.linalg.lstsq(A, P[m, 2], rcond=None)
            r = P[:, 2] - np.c_[np.ones(len(P)), P[:, 0], P[:, 1]] @ coef
            s = 1.4826 * np.median(np.abs(r[m]))
            m = np.abs(r) < 4 * s + 0.05
        return coef

    def footprint(self, name, plane, step=64):
        """Ground footprint polygon (n,2) of an image on the plane."""
        c = self.cams[name]
        u = np.r_[np.arange(0, c.w, step), np.full(c.h // step, c.w - 1),
                  np.arange(c.w - 1, 0, -step), np.zeros(c.h // step)]
        v = np.r_[np.zeros(len(np.arange(0, c.w, step))), np.arange(0, c.h, step)[:c.h // step],
                  np.full(len(np.arange(c.w - 1, 0, -step)), c.h - 1),
                  np.arange(c.h - 1, 0, -step)[:c.h // step]]
        return ray_plane(c, np.c_[u, v], plane)[:, :2]


def ray_plane(cam, uv, plane):
    """Intersect pixel rays with z = a + b x + c y. uv (n,2) -> (n,3)."""
    a, b, cc = plane
    d = np.c_[(uv[:, 0] - cam.cx) / cam.f, (uv[:, 1] - cam.cy) / cam.f, np.ones(len(uv))] @ cam.R
    # C + s d on plane: Cz + s dz = a + b (Cx + s dx) + c (Cy + s dy)
    s = (a + b * cam.C[0] + cc * cam.C[1] - cam.C[2]) / (d[:, 2] - b * d[:, 0] - cc * d[:, 1])
    return cam.C + s[:, None] * d


def rect_rotation(c1, c2):
    """Common rectifying rotation (world -> rect): x along the baseline,
    z close to the mean viewing direction."""
    x = c2.C - c1.C
    x /= np.linalg.norm(x)
    z0 = c1.R[2] + c2.R[2]                  # optical axes (world) of both
    y = np.cross(z0, x); y /= np.linalg.norm(y)
    z = np.cross(x, y)
    return np.vstack([x, y, z])
