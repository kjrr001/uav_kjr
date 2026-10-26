#!/usr/bin/env python3
"""ODM-KJR without ODM, step 0: initial cameras from the DJI EXIF / XMP tags.

  position  GPSLatitude / GPSLongitude -> UTM (zone from the longitude), AbsoluteAltitude
  rotation  GimbalYaw/Pitch/RollDegree (DJI; yaw clockwise from north, pitch -90 = nadir);
            else FlightYawDegree; else the direction of travel along the GPS track (image top
            forward), nadir
  intrinsics focal, first available: CalibratedFocalLength (px, DJI factory) ->
            FocalLengthIn35mmFormat (35 mm diagonal 43.27 mm) -> FocalLength x
            FocalPlaneXResolution; principal point CalibratedOpticalCenterX,Y else the centre;
            brown k1..k3, p1, p2 = 0 (self-calibrated later)
  heights   AbsoluteAltitude (DJI) else GPSAltitude; ground = AbsoluteAltitude -
            RelativeAltitude, else --ground-z, else flying height --agl
Writes an ODM-like run dir that every kjr step reads:
  opensfm/reconstruction.json            (brown camera, shots, ground points for the plane)
  opensfm/undistorted/reconstruction.json (perspective, same shots)
  odm_georeferencing/coords.txt          (UTM zone, integer offset)
  odm_orthophoto/odm_orthophoto.tif      (empty sparse grid: extent of the footprints, GSD)
  kjr_exif.json                          (raw tags, EPSG)
World frame = UTM minus the offset (as ODM / OpenSfM). Heights stay the GPS ones until the
GCP alignment of bba.py --align-gcp.
Usage: exif_init.py <image dir> <out dir> [--ground-z Z]
"""
import argparse, json, os, subprocess

import numpy as np

TAGS = ["GPSLatitude", "GPSLongitude", "GPSAltitude", "AbsoluteAltitude", "RelativeAltitude",
        "GimbalYawDegree", "GimbalPitchDegree", "GimbalRollDegree", "FlightYawDegree",
        "CalibratedFocalLength", "CalibratedOpticalCenterX", "CalibratedOpticalCenterY",
        "FocalLength", "FocalLengthIn35mmFormat", "FocalPlaneXResolution", "FocalPlaneResolutionUnit",
        "ImageWidth", "ImageHeight", "Model", "DateTimeOriginal"]


def focal_px(e, w, h):
    """focal in px and where it came from"""
    if e.get("CalibratedFocalLength"):
        return float(e["CalibratedFocalLength"]), "CalibratedFocalLength"
    if e.get("FocalLengthIn35mmFormat"):
        return float(e["FocalLengthIn35mmFormat"]) * np.hypot(w, h) / 43.267, "FocalLengthIn35mmFormat"
    if e.get("FocalLength") and e.get("FocalPlaneXResolution"):
        unit = {2: 25.4, 3: 10.0, 4: 1.0, 5: 0.001}.get(int(e.get("FocalPlaneResolutionUnit", 2)), 25.4)
        return float(e["FocalLength"]) * float(e["FocalPlaneXResolution"]) / unit, "FocalLength x FocalPlaneXResolution"
    raise SystemExit("no focal length in EXIF (CalibratedFocalLength / FocalLengthIn35mmFormat / FocalLength)")


def yaws(ex, E):
    """per image yaw (deg, clockwise from north) and the source"""
    if all("GimbalYawDegree" in e for e in ex):
        return np.array([float(e["GimbalYawDegree"]) for e in ex]), "GimbalYawDegree"
    if all("FlightYawDegree" in e for e in ex):
        return np.array([float(e["FlightYawDegree"]) for e in ex]), "FlightYawDegree"
    order = np.argsort([e.get("DateTimeOriginal", "") for e in ex], kind="stable")
    y = np.zeros(len(ex))
    P = E[order]
    for k in range(len(P)):
        a, b = max(k - 1, 0), min(k + 1, len(P) - 1)
        d = P[b] - P[a]
        y[order[k]] = np.degrees(np.arctan2(d[0], d[1]))
    return y, "GPS track (image top = direction of travel)"


def read_exif(img_dir):
    out = subprocess.run(["exiftool", "-n", "-j", "-q", "-ext", "JPG", "-ext", "jpg"]
                         + [f"-{t}" for t in TAGS] + [img_dir], capture_output=True, text=True, check=True)
    return sorted(json.loads(out.stdout), key=lambda d: d["SourceFile"])


def load_camera(path, w, h):
    """brown intrinsics (OpenSfM normalised units) of the camera with this image size"""
    d = json.load(open(path))
    if isinstance(d, list):                      # reconstruction.json
        d = d[0]["cameras"]
    if "intrinsics" in d:                        # kjr bba.json
        d = {"bba": d["intrinsics"]}
    for c in d.values():
        if int(c.get("width", w)) == w and int(c.get("height", h)) == h:
            keys = ["focal_x", "c_x", "c_y", "k1", "k2", "k3", "p1", "p2"]
            out = {k: float(c.get(k, c.get("focal", 0.0) if k == "focal_x" else 0.0)) for k in keys}
            out["focal_y"] = out["focal_x"]
            return out
    raise SystemExit(f"{path}: no camera of {w}x{h}")


def rot_from_gimbal(yaw, pitch, roll):
    """world (E, N, U) -> camera (x right, y down, z forward), OpenSfM convention.
    Nadir (pitch -90) with yaw psi: image top = heading (sin psi, cos psi, 0)."""
    psi = np.radians(yaw)
    h = np.array([np.sin(psi), np.cos(psi), 0.0])          # heading
    r = np.array([np.cos(psi), -np.sin(psi), 0.0])          # right of the heading
    # pitch: tilt from nadir about the right axis (pitch -90 -> optical axis down)
    tilt = np.radians(pitch + 90.0)
    z = np.array([0.0, 0.0, -1.0]) * np.cos(tilt) + h * np.sin(tilt)
    y = np.cross(z, r)
    x = r
    R = np.vstack([x, y, z])
    if roll:                                                  # about the optical axis
        c, s = np.cos(np.radians(roll)), np.sin(np.radians(roll))
        R = np.array([[c, s, 0], [-s, c, 0], [0, 0, 1]]) @ R
    return R


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("images"); ap.add_argument("out")
    ap.add_argument("--ground-z", type=float, default=None,
                    help="ground height for footprints (default: AbsoluteAltitude - RelativeAltitude)")
    ap.add_argument("--agl", type=float, default=None, help="flying height above ground, m (no RelativeAltitude)")
    ap.add_argument("--camera", default=None,
                    help="camera calibration to start from (and hold with bba.py --intr-prior): an"
                         " OpenSfM/ODM cameras.json or reconstruction.json, or a kjr bba.json")
    a = ap.parse_args()
    import cv2
    from osgeo import osr, gdal
    ex = read_exif(a.images)
    names = [os.path.basename(e["SourceFile"]) for e in ex]
    lat = np.array([e["GPSLatitude"] for e in ex]); lon = np.array([e["GPSLongitude"] for e in ex])
    zone = int((lon.mean() + 180) // 6) + 1
    epsg = (32600 if lat.mean() >= 0 else 32700) + zone
    src = osr.SpatialReference(); src.ImportFromEPSG(4326)
    dst = osr.SpatialReference(); dst.ImportFromEPSG(epsg)
    for s in (src, dst):
        s.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    tr = osr.CoordinateTransformation(src, dst)
    E = np.array([tr.TransformPoint(float(x), float(y))[:2] for x, y in zip(lon, lat)])
    off = np.round(E.mean(0))
    alt = np.array([e.get("AbsoluteAltitude", e.get("GPSAltitude", np.nan)) for e in ex], float)
    rel = np.array([e.get("RelativeAltitude", np.nan) for e in ex], float)
    if a.ground_z is not None:
        zg = a.ground_z
    elif np.isfinite(rel).all():
        zg = float(np.nanmedian(alt - rel))
    elif a.agl is not None:
        zg = float(np.median(alt)) - a.agl
    else:
        raise SystemExit("no RelativeAltitude: give --ground-z or --agl")
    C = np.c_[E - off, alt]
    w, h = int(ex[0]["ImageWidth"]), int(ex[0]["ImageHeight"])
    S = max(w, h)
    if len({(e.get("Model"), e["ImageWidth"], e["ImageHeight"], e.get("FocalLength")) for e in ex}) > 1:
        raise SystemExit("more than one camera / focal in the block: not supported (one camera model)")
    fpx, fsrc = focal_px(ex[0], w, h)
    cxp = float(ex[0].get("CalibratedOpticalCenterX", w / 2)); cyp = float(ex[0].get("CalibratedOpticalCenterY", h / 2))
    yaw, ysrc = yaws(ex, E)
    model = ex[0].get("Model", "camera")
    key = f"v2 {model.lower()} {w} {h} brown exif"
    cam = {"projection_type": "brown", "width": w, "height": h, "focal_x": fpx / S, "focal_y": fpx / S,
           "c_x": (cxp - w / 2) / S, "c_y": (cyp - h / 2) / S, "k1": 0.0, "k2": 0.0, "k3": 0.0, "p1": 0.0, "p2": 0.0}
    if a.camera:
        cam.update(load_camera(a.camera, w, h))
        fpx = cam["focal_x"] * S
        fsrc = f"--camera {a.camera}"
    shots, ushots = {}, {}
    for k, (n, e) in enumerate(zip(names, ex)):
        R = rot_from_gimbal(yaw[k], float(e.get("GimbalPitchDegree", -90.0)), float(e.get("GimbalRollDegree", 0) or 0))
        rv = cv2.Rodrigues(R)[0].ravel()
        t = -R @ C[k]
        s = {"rotation": rv.tolist(), "translation": t.tolist(), "camera": key,
             "gps_position": C[k].tolist(), "gps_dop": 5.0}
        shots[n] = s; ushots[n + ".tif"] = dict(s)
    # footprints on z = zg -> ground points (plane fit) and the grid extent
    pts, allxy = {}, []
    for k, n in enumerate(names):
        R = cv2.Rodrigues(np.array(shots[n]["rotation"]))[0]
        uv = np.array([[0, 0], [w, 0], [w, h], [0, h], [w / 2, h / 2]], float)
        d = np.c_[(uv - [cxp, cyp]) / fpx, np.ones(len(uv))] @ R
        X = C[k] + ((zg - C[k, 2]) / d[:, 2])[:, None] * d
        allxy.append(X[:4, :2])
        pts[str(k)] = {"coordinates": X[4].tolist(), "color": [0, 0, 0]}
    rec = [{"cameras": {key: cam}, "shots": shots, "points": pts,
            "reference_lla": {"latitude": float(lat.mean()), "longitude": float(lon.mean()), "altitude": 0.0}}]
    ucam = {"projection_type": "perspective", "width": w, "height": h, "focal": fpx / S, "k1": 0.0, "k2": 0.0}
    urec = [{"cameras": {key: ucam}, "shots": ushots, "points": pts}]
    os.makedirs(os.path.join(a.out, "opensfm/undistorted"), exist_ok=True)
    os.makedirs(os.path.join(a.out, "odm_georeferencing"), exist_ok=True)
    os.makedirs(os.path.join(a.out, "odm_orthophoto"), exist_ok=True)
    json.dump(rec, open(os.path.join(a.out, "opensfm/reconstruction.json"), "w"))
    json.dump(urec, open(os.path.join(a.out, "opensfm/undistorted/reconstruction.json"), "w"))
    hemi = "N" if lat.mean() >= 0 else "S"
    open(os.path.join(a.out, "odm_georeferencing/coords.txt"), "w").write(
        f"WGS84 UTM {zone}{hemi}\n{int(off[0])} {int(off[1])}\n")
    # empty grid: GSD of the median flying height, extent of all footprints (+20 m)
    gsd = float(np.median(C[:, 2] - zg)) / fpx
    xy = np.concatenate(allxy) + off
    x0, y0 = np.floor(xy.min(0) - 20); x1, y1 = np.ceil(xy.max(0) + 20)
    W, H = int((x1 - x0) / gsd), int((y1 - y0) / gsd)
    d = gdal.GetDriverByName("GTiff").Create(os.path.join(a.out, "odm_orthophoto/odm_orthophoto.tif"), W, H, 1,
                                             gdal.GDT_Byte, ["SPARSE_OK=TRUE", "TILED=YES"])
    d.SetGeoTransform((x0, gsd, 0, y1, 0, -gsd)); d.SetProjection(dst.ExportToWkt()); d = None
    json.dump({"epsg": epsg, "offset": off.tolist(), "ground_z": zg, "gsd": gsd, "exif": ex},
              open(os.path.join(a.out, "kjr_exif.json"), "w"), indent=0)
    print(f"{len(names)} images, {model} {w}x{h} f {fpx:.1f} px ({fsrc}), yaw from {ysrc}; EPSG {epsg}, offset {off}; "
          f"ground z {zg:.1f}, flying height {np.median(C[:, 2] - zg):.1f} m, GSD {gsd*100:.2f} cm, grid {W}x{H}")


if __name__ == "__main__":
    main()
