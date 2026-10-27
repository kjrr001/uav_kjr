#!/usr/bin/env python3
"""ODM-KJR without ODM: import an upstream OpenSfM (mapillary/OpenSfM, BSD) reconstruction.

OpenSfM is run on the images and the GCP list by the user (tools/run_opensfm.sh); it leaves
reconstruction.json in a topocentric (ENU) frame about reference_lla. This writes the
ODM-like run dir every kjr step reads (as exif_init.py does):
  opensfm/reconstruction.json             brown camera, shots, points in UTM - offset
  opensfm/undistorted/reconstruction.json perspective camera (focal = brown focal), same shots
  odm_georeferencing/coords.txt           UTM zone, integer offset
  odm_orthophoto/odm_orthophoto.tif       empty sparse grid (extent of the points, GSD)
  odm_dem/dsm.tif                         reference DSM from the OpenSfM points (refdsm.py)
  kjr_exif.json                           EPSG, offset
ENU -> UTM: points and camera centres exact (ENU -> ECEF -> UTM; heights in the GCP datum,
which OpenSfM took as ellipsoidal both ways); camera rotations from a similarity fitted per shot
to the exact transform near it (one similarity for a 2 km block is off by ~0.1 m: curvature).
Points outside the 0.5-99.5 % box (+50 m) are left out of the grid extent and reference DSM.
Only the largest reconstruction is used; the shots of the others are listed.
Usage: osfm_import.py <opensfm dataset dir> <out run dir> [--cell 2] [--gsd M]
"""
import argparse, json, os, subprocess, sys

import cv2
import numpy as np
from osgeo import gdal, osr

gdal.UseExceptions()
A_WGS, E2 = 6378137.0, 6.69437999014e-3


def enu_to_ecef(P, lat, lon, alt):
    la, lo = np.radians(lat), np.radians(lon)
    N = A_WGS / np.sqrt(1 - E2 * np.sin(la) ** 2)
    o = np.array([(N + alt) * np.cos(la) * np.cos(lo), (N + alt) * np.cos(la) * np.sin(lo),
                  (N * (1 - E2) + alt) * np.sin(la)])
    R = np.array([[-np.sin(lo), -np.sin(la) * np.cos(lo), np.cos(la) * np.cos(lo)],
                  [np.cos(lo), -np.sin(la) * np.sin(lo), np.cos(la) * np.sin(lo)],
                  [0, np.cos(la), np.sin(la)]])
    return P @ R.T + o


def umeyama(X, Y):
    """Y ~ s Q X + T"""
    mx, my = X.mean(0), Y.mean(0)
    Xc, Yc = X - mx, Y - my
    U, S, Vt = np.linalg.svd(Yc.T @ Xc / len(X))
    D = np.diag([1, 1, np.sign(np.linalg.det(U @ Vt))])
    Q = U @ D @ Vt
    s = np.trace(np.diag(S) @ D) / Xc.var(0).sum()
    return s, Q, my - s * Q @ mx


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("osfm"); ap.add_argument("out")
    ap.add_argument("--cell", type=float, default=2.0, help="reference DSM cell (m)")
    ap.add_argument("--gsd", type=float, default=None, help="ortho grid cell (m); default flying height / focal")
    a = ap.parse_args()
    recs = json.load(open(os.path.join(a.osfm, "reconstruction.json")))
    recs.sort(key=lambda r: -len(r["shots"]))
    rec = recs[0]
    for r in recs[1:]:
        print(f"WARNING: second reconstruction with {len(r['shots'])} shots not used: {sorted(r['shots'])[:5]} ...")
    nimg = len([f for f in os.listdir(os.path.join(a.osfm, "images"))])
    print(f"reconstruction: {len(rec['shots'])} of {nimg} images, {len(rec['points'])} points")
    if len(rec["cameras"]) != 1:
        raise SystemExit("more than one camera: not supported")
    ref = rec["reference_lla"]
    lat, lon, alt = ref["latitude"], ref["longitude"], ref["altitude"]
    zone = int((lon + 180) // 6) + 1
    hemi = "N" if lat >= 0 else "S"
    epsg = (32600 if lat >= 0 else 32700) + zone
    utm = osr.SpatialReference(); utm.ImportFromEPSG(epsg)
    ecef = osr.SpatialReference(); ecef.ImportFromEPSG(4978)
    for s in (utm, ecef):
        s.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    tr = osr.CoordinateTransformation(ecef, utm)

    def to_utm(P):
        return np.array(tr.TransformPoints(enu_to_ecef(np.asarray(P, float).reshape(-1, 3), lat, lon, alt).tolist()))

    names = sorted(rec["shots"])
    Rs, Cs = [], []
    for n in names:
        s = rec["shots"][n]
        R = cv2.Rodrigues(np.array(s["rotation"], float))[0]
        Rs.append(R); Cs.append(-R.T @ np.array(s["translation"], float))
    Cs = np.array(Cs)
    pts = np.array([p["coordinates"] for p in rec["points"].values()], float)
    # per shot: centre exactly, rotation from a similarity fitted to the exact transform over the
    # column from the camera down to the ground (one similarity for the block is off by ~0.1 m
    # at 2 km from the earth curvature)
    off = np.round(to_utm(Cs)[:, :2].mean(0))
    o3 = np.array([off[0], off[1], 0.0])
    zg0 = float(np.median(pts[:, 2]))
    u = np.linspace(-1, 1, 3)
    shots, worst = {}, 0.0
    for n, R, C in zip(names, Rs, Cs):
        hgt = max(C[2] - zg0, 1.0)
        g = C + np.stack(np.meshgrid(u * hgt, u * hgt, np.linspace(-hgt, 0, 3), indexing="ij"), -1).reshape(-1, 3)
        G = to_utm(g) - o3
        s_, Q, T = umeyama(g, G)
        worst = max(worst, np.abs(s_ * g @ Q.T + T - G).max())
        R2 = R @ Q.T
        C2 = to_utm(C)[0] - o3
        d = dict(rec["shots"][n])
        d["rotation"] = cv2.Rodrigues(R2)[0].ravel().tolist()
        d["translation"] = (-R2 @ C2).tolist()
        if "gps_position" in d:
            d["gps_position"] = (to_utm(d["gps_position"])[0] - o3).tolist()
        shots[n] = d
    print(f"ENU -> UTM {epsg}: centres and points exact, rotations from local similarities "
          f"(max residual {worst * 1000:.2f} mm)")
    P2 = to_utm(pts) - o3
    points = {k: {"coordinates": P2[i].tolist(), "color": p.get("color", [0, 0, 0])}
              for i, (k, p) in enumerate(rec["points"].items())}
    key, cam = next(iter(rec["cameras"].items()))
    w, h = cam["width"], cam["height"]
    S = max(w, h)
    f = cam.get("focal_x", cam.get("focal")) if cam["projection_type"] == "brown" else cam["focal"]
    fy = cam.get("focal_y", f)
    ucam = {"projection_type": "perspective", "width": w, "height": h, "focal": (f + fy) / 2, "k1": 0.0, "k2": 0.0}
    out = {"cameras": {key: cam}, "shots": shots, "points": points, "reference_lla": ref}
    uout = {"cameras": {key: ucam}, "shots": {n + ".tif": shots[n] for n in names}, "points": points}
    for d in ("opensfm/undistorted", "odm_georeferencing", "odm_orthophoto", "odm_dem"):
        os.makedirs(os.path.join(a.out, d), exist_ok=True)
    json.dump([out], open(os.path.join(a.out, "opensfm/reconstruction.json"), "w"))
    json.dump([uout], open(os.path.join(a.out, "opensfm/undistorted/reconstruction.json"), "w"))
    open(os.path.join(a.out, "odm_georeferencing/coords.txt"), "w").write(f"WGS84 UTM {zone}{hemi}\n{int(off[0])} {int(off[1])}\n")
    # sparse ortho grid: points extent + 20 m, GSD at the median flying height
    q = np.percentile(P2, [0.5, 99.5], axis=0)
    inb = np.all((P2 >= q[0] - [50, 50, 20]) & (P2 <= q[1] + [50, 50, 20]), 1)
    zg = np.median(P2[inb, 2])
    C2 = np.array([-cv2.Rodrigues(np.array(shots[n]["rotation"]))[0].T @ np.array(shots[n]["translation"]) for n in names])
    gsd = a.gsd or float(np.median(C2[:, 2]) - zg) / (f * S)
    # extent: image footprints on the median ground plane (+20 m), not only the tie points
    # (the 0.5-99.5 % point box cut the block edges)
    uv = (np.array([[0, 0], [w, 0], [w, h], [0, h]], float) - [w / 2, h / 2]) / (f * S)
    xy = []
    for n, C in zip(names, C2):
        R = cv2.Rodrigues(np.array(shots[n]["rotation"]))[0]
        dv = np.c_[uv, np.ones(4)] @ R
        xy.append((C + ((zg - C[2]) / dv[:, 2])[:, None] * dv)[:, :2])
    xy = np.concatenate(xy) + off
    x0, y0 = np.floor(xy.min(0) - 20); x1, y1 = np.ceil(xy.max(0) + 20)
    W, H = int((x1 - x0) / gsd), int((y1 - y0) / gsd)
    d = gdal.GetDriverByName("GTiff").Create(os.path.join(a.out, "odm_orthophoto/odm_orthophoto.tif"), W, H, 1,
                                             gdal.GDT_Byte, ["SPARSE_OK=TRUE", "TILED=YES"])
    d.SetGeoTransform((x0, gsd, 0, y1, 0, -gsd)); d.SetProjection(utm.ExportToWkt()); d = None
    json.dump({"epsg": epsg, "offset": off.tolist(), "ground_z": float(zg), "gsd": gsd, "source": "opensfm"},
              open(os.path.join(a.out, "kjr_exif.json"), "w"), indent=0)
    print(f"{len(names)} shots, f {f * S:.1f} px, EPSG {epsg}, offset {off}, GSD {gsd * 100:.2f} cm, grid {W}x{H}")
    # reference DSM from the OpenSfM points
    pd = os.path.join(a.out, "osfm_points"); os.makedirs(pd, exist_ok=True)
    np.save(os.path.join(pd, "points.npy"), P2[inb].astype(np.float32))
    subprocess.check_call([sys.executable, os.path.join(os.path.dirname(os.path.abspath(__file__)), "refdsm.py"),
                           pd, a.out, "--cell", str(a.cell)])


if __name__ == "__main__":
    main()
