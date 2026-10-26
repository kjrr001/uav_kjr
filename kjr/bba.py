#!/usr/bin/env python3
"""ODM-KJR BBA step 2: bundle adjustment with self-calibration, from the
sub-pixel ties of kjr/ties.py, starting from the ODM/OpenSfM solution.

Why: the KJ-1 DSM ripple (rings, 1-3 cm) is fixed to IMAGE coordinates and is
the same in ODM's own DSM (OpenMVS) and in ours - it comes from the shared
camera model, not from any matcher (diag_imgfixed.py, 2026-10-25). And the
east strip has 10-17 cm steps between consecutive pairs (pose error).

Camera model = OpenSfM brown (normalised by size = max(w, h)) + corrections in
pixels at the brown pixel position: a radial one, piecewise linear in r (node
every RSTEP px), and a 2-D bilinear grid (node every GSTEP px) for x and y. Both
are kept orthogonal to the smooth modes (radial: 1, r, r^3.., grid: degree <= 3),
which trade with block curvature:
  xn = Xc/Zc, r2 = xn^2 + yn^2, radial = 1 + k1 r2 + k2 r2^2 + k3 r2^3
  xd = xn radial + 2 p1 xn yn + p2 (r2 + 2 xn^2), yd likewise
  u0 = (f xd + cx) size + w/2 - 0.5 ; v0 likewise
  u = u0 + dr(rho) (u0 - w/2 + 0.5) / rho + gx(u0, v0), rho = |(u0, v0) - centre|  (v: gy)
World = ODM's (UTM - coords.txt offset), Xc = R(rotvec) X + t (OpenSfM shot).

Observations
  ties     multi-view tracks: the first image's SIFT position (unmoved by ALSC)
           groups all its pairs -> one 3-D point per anchor
  GCPs     gcp_list (original image px), fixed ground coordinates
  priors   poses to the ODM solution (weak), dr to 0 and smooth
Solver: Levenberg-Marquardt, points eliminated by the Schur complement,
Cauchy-weighted (IRLS) residuals.

Writes <out>/bba.json (intrinsics, dr nodes, shots, stats) and an ODM-like
run dir <out>/odm (opensfm/reconstruction.json, opensfm/undistorted/
reconstruction.json, odm_georeferencing/coords.txt) for kjr/undistort.py and
kjr/stereo.py.
Usage: bba.py <odm run> <ties dir> <gcp list> <out>
"""
import argparse
import glob
import json
import os
import shutil
import time

import numpy as np
from scipy import sparse
from scipy.spatial.transform import Rotation

RSTEP = 30.0           # px between dr nodes (the dr oscillation has a ~180-300 px period)
GSTEP = 64.0           # px between 2-D grid nodes (30 px would need a 16 GB dense Schur matrix)
SIG_G = 0.5            # px, grid prior to 0
SIG_G2 = 0.2           # px, grid second difference (x and y)
SIG_TIE = 0.1          # px, nominal tie sigma (ALSC)
SIG_GCP = 0.3          # px, hand marks. 0.7 let 2.3 M ties tilt KJ-1 (DSM at GCP 1 +16 cm)
SIG_ROT = 0.02         # rad, pose prior to ODM
SIG_T = 1.0            # m
SIG_DR = 1.0           # px, dr prior to 0
SIG_DR2 = 0.1          # px, dr second difference
SIG_FIX = 1e-9         # intrinsics held at the ODM values (unless --free-intr)
MAX_OBS = 4_000_000    # tie observations kept; JJ-3 (FC9313, 17 M) ran out of 30 GB. KJ/JJ1/JJ2 are below it
SIG_ORTH = 1e-4        # px: dr is kept orthogonal to 1, r, r^3, r^5, r^7 - the smooth radial
                       # part belongs to f, k1..k3 (otherwise they trade, dr -> several px)
INTR = ["focal_x", "c_x", "c_y", "k1", "k2", "k3", "p1", "p2"]


class Model:
    def __init__(self, w, h):
        self.w, self.h, self.S = w, h, float(max(w, h))
        self.c0 = np.array([w / 2 - 0.5, h / 2 - 0.5])
        self.nr = int(np.ceil(np.hypot(w, h) / 2 / RSTEP)) + 2
        self.gx = int(np.ceil((w - 1) / GSTEP)) + 1
        self.gy = int(np.ceil((h - 1) / GSTEP)) + 1
        self.ng = self.gx * self.gy
        self.nq = self.nr + 2 * self.ng      # correction vector q = [dr, grid x, grid y]

    def brown(self, Xc, ip):
        f, cx, cy, k1, k2, k3, p1, p2 = ip
        xn, yn = Xc[:, 0] / Xc[:, 2], Xc[:, 1] / Xc[:, 2]
        r2 = xn * xn + yn * yn
        rad = 1 + k1 * r2 + k2 * r2 ** 2 + k3 * r2 ** 3
        xd = xn * rad + 2 * p1 * xn * yn + p2 * (r2 + 2 * xn * xn)
        yd = yn * rad + p1 * (r2 + 2 * yn * yn) + 2 * p2 * xn * yn
        return np.c_[(f * xd + cx) * self.S, (f * yd + cy) * self.S] + self.c0

    def dr_terms(self, u0):
        d = u0 - self.c0
        rho = np.maximum(np.hypot(d[:, 0], d[:, 1]), 1e-6)
        x = rho / RSTEP
        i0 = np.minimum(x.astype(int), self.nr - 2)
        a = x - i0
        return d / rho[:, None], i0, a

    def corr_terms(self, u0):
        """correction = sum_k q[idx[:, k]] * val[:, :, k]  (n, 2) px, linear in q.
        k: 2 radial nodes (along e), 4 grid nodes for x, 4 for y (bilinear)."""
        n = len(u0)
        e, i0, a = self.dr_terms(u0)
        gx = np.clip(u0[:, 0] / GSTEP, 0, self.gx - 1 - 1e-9); gy = np.clip(u0[:, 1] / GSTEP, 0, self.gy - 1 - 1e-9)
        jx, jy = gx.astype(int), gy.astype(int); fx, fy = gx - jx, gy - jy
        node = np.stack([jy * self.gx + jx, jy * self.gx + jx + 1, (jy + 1) * self.gx + jx, (jy + 1) * self.gx + jx + 1], 1)
        wb = np.stack([(1 - fx) * (1 - fy), fx * (1 - fy), (1 - fx) * fy, fx * fy], 1)
        idx = np.concatenate([np.c_[i0, i0 + 1], self.nr + node, self.nr + self.ng + node], 1)
        val = np.zeros((n, 2, 10))
        val[:, :, 0] = (1 - a)[:, None] * e; val[:, :, 1] = a[:, None] * e
        val[:, 0, 2:6] = wb; val[:, 1, 6:10] = wb
        return idx, val

    def corr(self, u0, q):
        idx, val = self.corr_terms(u0)
        return np.einsum("nck,nk->nc", val, q[idx])

    def project(self, Xc, ip, q):
        u0 = self.brown(Xc, ip)
        return u0 + self.corr(u0, q)


def load_schur():
    """C++ point-elimination kernel (kjr/csrc/schur.cpp); None -> scipy sparse products."""
    import ctypes
    so = os.path.join(os.path.dirname(os.path.abspath(__file__)), "_schur.so")
    if not os.path.exists(so):
        src = os.path.join(os.path.dirname(so), "csrc", "schur.cpp")
        if os.system(f"g++ -O3 -march=native -fopenmp -shared -fPIC {src} -o {so}") != 0:
            return None
    lib = ctypes.CDLL(so)
    P = np.ctypeslib.ndpointer
    lib.schur_points.argtypes = [ctypes.c_int64, P(np.int64), ctypes.c_int, P(np.int64), P(np.float64),
                                 P(np.float64), P(np.float64), ctypes.c_int64, P(np.float64), P(np.float64),
                                 P(np.float64), P(np.float64), ctypes.c_int]
    lib.schur_points.restype = ctypes.c_int
    return lib


def rot(rv):
    return Rotation.from_rotvec(rv).as_matrix()


def cam_xc(R, t, cam, X):
    return np.einsum("nij,nj->ni", R[cam], X) + t[cam]


def undistort_px(model, ip, q, uv, it=20):
    """pixel -> normalised (xn, yn), fixed point on the brown model (+ correction q)."""
    f, cx, cy = ip[:3]
    u0 = uv.copy()
    for _ in range(3):                       # correction is evaluated at the brown position
        u0 = uv - model.corr(u0, q)
    uv = u0
    xd = ((uv - model.c0) / model.S - [cx, cy]) / f
    xn = xd.copy()
    for _ in range(it):
        p = model.brown(np.c_[xn, np.ones(len(xn))], ip)
        xn -= ((p - model.c0) / model.S - [cx, cy]) / f - xd
    return xn


def triangulate(R, t, model, ip, dr, cam, uv, pt, npt):
    xn = undistort_px(model, ip, dr, uv)
    d = np.einsum("nji,nj->ni", R[cam], np.c_[xn, np.ones(len(xn))])   # R^T [xn 1]
    d /= np.linalg.norm(d, axis=1)[:, None]
    C = -np.einsum("nji,nj->ni", R[cam], t[cam])
    M = np.eye(3)[None] - d[:, :, None] * d[:, None, :]
    A = np.zeros((npt, 3, 3)); b = np.zeros((npt, 3))
    np.add.at(A, pt, M); np.add.at(b, pt, np.einsum("nij,nj->ni", M, C))
    A += 1e-9 * np.eye(3)
    return np.linalg.solve(A, b[..., None])[..., 0]


def load_ties(tie_dir, names, max_obs=MAX_OBS):
    """multi-view tracks: anchor = (first image, exact SIFT position).
    Above max_obs, anchors are thinned by a hash of their position: the same anchor is kept
    or dropped in every pair, so multi-view tracks stay whole."""
    idx = {n: k for k, n in enumerate(names)}
    files = [f for f in sorted(glob.glob(os.path.join(tie_dir, "*.npz")))]
    nobs = 0
    for fn in files:
        with np.load(fn) as d:
            n = len(d["pa"])
        nobs += 2 * n if n >= 30 else 0
    frac = min(1.0, max_obs / max(nobs, 1))
    obs_cam, obs_uv, obs_pt = [], [], []
    anchors = {}
    npt = 0
    for fn in files:
        a, b = os.path.basename(fn)[:-4].split("__")
        with np.load(fn) as d:
            pa, pb = d["pa"], d["pb"]
        if len(pa) < 30:
            continue
        if frac < 1:
            h = np.sin(pa[:, 0] * 12.9898 + pa[:, 1] * 78.233) * 43758.5453
            sel = (h - np.floor(h)) < frac
            pa, pb = pa[sel], pb[sel]
        ia, ib = idx[a], idx[b]
        tab = anchors.setdefault(ia, {})
        for p, q in zip(map(tuple, pa), pb):
            k = tab.get(p)
            if k is None:
                k = tab[p] = npt; npt += 1
                obs_cam.append(ia); obs_uv.append(p); obs_pt.append(k)
            obs_cam.append(ib); obs_uv.append(q); obs_pt.append(k)
    cam, uv, pt = np.array(obs_cam), np.array(obs_uv, float), np.array(obs_pt)
    # a point seen twice in one image (two anchors merged wrongly) cannot happen
    # with star tracks; drop duplicates (same point, same image) anyway
    key = pt * 1000 + cam
    _, first = np.unique(key, return_index=True)
    keep = np.zeros(len(pt), bool); keep[first] = True
    return cam[keep], uv[keep], pt[keep], npt


def load_tracks(path, names, max_obs=MAX_OBS):
    """multi-view tracks of kjr/tracks.py (one ground point per track, every observation
    measured from the track's reference). Above max_obs whole tracks are thinned."""
    d = np.load(path)
    idx = {n: k for k, n in enumerate(names)}
    tn = [str(n) for n in d["names"]]
    m = np.array([idx.get(n, -1) for n in tn])
    cam, uv, pt = m[d["cam"]], d["uv"], d["pt"]
    keep = cam >= 0
    if keep.sum() > max_obs:
        h = np.sin(pt * 12.9898 + 78.233) * 43758.5453
        keep &= (h - np.floor(h)) < max_obs / keep.sum()
    cam, uv, pt = cam[keep], uv[keep], pt[keep]
    ok = np.bincount(pt)[pt] >= 2
    cam, uv, pt = cam[ok], uv[ok], pt[ok]
    _, pt = np.unique(pt, return_inverse=True)
    return cam, uv, pt, pt.max() + 1


def align_to_gcp(rv, t, model, ip, g_cam, g_uv, g_X, g_id, say):
    """Similarity (scale, rotation, shift; Umeyama) taking the GCPs triangulated from the
    input cameras onto the surveyed coordinates, applied to the cameras."""
    ug = sorted(set(g_id)); gi = np.array([ug.index(x) for x in g_id])
    nv = np.bincount(gi, minlength=len(ug))
    R = rot(rv)
    Xg = triangulate(R, t, model, ip, np.zeros(model.nq), g_cam, g_uv, gi, len(ug))
    Xs = np.array([g_X[gi == k][0] for k in range(len(ug))])
    ok = nv >= 2
    if ok.sum() < 3:
        say(f"align-gcp: only {ok.sum()} GCPs seen in 2+ images - not aligned")
        return rv, t
    A, B = Xg[ok], Xs[ok]
    ma, mb = A.mean(0), B.mean(0)
    U, S, Vt = np.linalg.svd((B - mb).T @ (A - ma))
    D = np.diag([1, 1, np.sign(np.linalg.det(U @ Vt))])
    Q = U @ D @ Vt
    sc = np.trace(np.diag(S) @ D) / ((A - ma) ** 2).sum()
    T = mb - sc * Q @ ma
    res = (sc * (Q @ A.T)).T + T - B
    C = -np.einsum("nji,nj->ni", R, t)
    C2 = sc * C @ Q.T + T
    R2 = R @ Q.T
    t2 = -np.einsum("nij,nj->ni", R2, C2)
    rv2 = Rotation.from_matrix(R2).as_rotvec()
    say(f"align-gcp: {ok.sum()} GCPs, scale {sc:.5f}, rotation {np.degrees(Rotation.from_matrix(Q).magnitude()):.3f} deg, "
        f"mean camera shift {np.round((C2 - C).mean(0), 2)} m; residual RMS {np.sqrt((res ** 2).sum(1).mean()):.2f} m")
    return rv2, t2


def main():
    global SIG_T, SIG_ROT
    ap = argparse.ArgumentParser()
    ap.add_argument("odm"); ap.add_argument("ties"); ap.add_argument("gcp"); ap.add_argument("out")
    ap.add_argument("--iters", type=int, default=12)
    ap.add_argument("--sig-t", type=float, default=SIG_T, help="m, camera centre prior to the input solution")
    ap.add_argument("--sig-rot", type=float, default=SIG_ROT, help="rad, rotation prior to the input solution")
    ap.add_argument("--no-gcp", action="store_true", help="ties only (GCP list ignored)")
    ap.add_argument("--align-gcp", action="store_true",
                    help="first move the input cameras by the similarity that maps the GCPs triangulated"
                         " from them onto the surveyed ones (input from GPS/EXIF: datum, heights)")
    ap.add_argument("--tracks", default=None,
                    help="multi-view tracks (kjr/tracks.py tracks.npz) instead of the pairwise ties")
    ap.add_argument("--no-dr", action="store_true", help="brown only (no radial correction)")
    ap.add_argument("--tol", type=float, default=1e-5, help="stop when the relative cost decrease is below this")
    ap.add_argument("--sparse-schur", action="store_true",
                    help="reference path: Schur complement with scipy sparse products (3-4x slower)")
    ap.add_argument("--free-intr", action="store_true",
                    help="also adjust f, cx, cy, k1-k3, p1, p2. Off by default: on KJ-1 (parallel nadir"
                         " strips, 5 GCPs) f and k1-k3 trade with block curvature -> a 1 m dome")
    args = ap.parse_args()
    SIG_T, SIG_ROT = args.sig_t, args.sig_rot
    os.makedirs(args.out, exist_ok=True)
    log = open(os.path.join(args.out, "bba.log"), "a")

    def say(*a):
        s = " ".join(str(x) for x in a); print(s, flush=True); log.write(s + "\n"); log.flush()

    rec = json.load(open(os.path.join(args.odm, "opensfm/reconstruction.json")))[0]
    cam_key, camd = next(iter(rec["cameras"].items()))
    names = sorted(rec["shots"])
    model = Model(camd["width"], camd["height"])
    ip0 = np.array([camd[k] for k in INTR], float)
    rv0 = np.array([rec["shots"][n]["rotation"] for n in names], float)
    t0 = np.array([rec["shots"][n]["translation"] for n in names], float)
    with open(os.path.join(args.odm, "odm_georeferencing/coords.txt")) as f:
        f.readline(); off = np.array(list(map(float, f.readline().split()[:2])))
    nc = len(names)

    # ---------------- observations
    cam, uv, pt, npt = load_tracks(args.tracks, names) if args.tracks else load_ties(args.ties, names)
    say(f"ties: {npt} points, {len(cam)} observations, "
        f"views/point mean {len(cam) / npt:.2f}, >=3 views {np.mean(np.bincount(pt) >= 3) * 100:.0f}%")
    idx = {n: k for k, n in enumerate(names)}
    g_cam, g_uv, g_X, g_id = [], [], [], []
    for line in ([] if args.no_gcp else open(args.gcp)):
        s = line.split()
        if len(s) < 6 or s[5] not in idx:
            continue
        g_cam.append(idx[s[5]]); g_uv.append([float(s[3]), float(s[4])])
        g_X.append([float(s[0]) - off[0], float(s[1]) - off[1], float(s[2])]); g_id.append(s[6] if len(s) > 6 else "")
    g_cam = np.array(g_cam, int); g_uv = np.array(g_uv, float).reshape(-1, 2); g_X = np.array(g_X, float).reshape(-1, 3)
    say(f"GCP observations: {len(g_cam)}")
    if args.align_gcp and len(g_cam):
        rv0, t0 = align_to_gcp(rv0, t0, model, ip0, g_cam, g_uv, g_X, g_id, say)

    dr = np.zeros(model.nq)                  # all corrections (name kept: dr = q)
    rv, t, ip = rv0.copy(), t0.copy(), ip0.copy()
    R = rot(rv)
    X = triangulate(R, t, model, ip, dr, cam, uv, pt, npt)
    # wrong matches with near-parallel rays triangulate behind a camera or km away; their
    # projection is NaN/inf and stops the adjustment (JJ-3: 351 of 17 M observations)
    zc = cam_xc(R, t, cam, X[pt])[:, 2]
    zm = np.median(zc)
    badpt = np.unique(pt[(zc < 0.5 * zm) | (zc > 2.0 * zm)])
    if len(badpt):
        keep = ~np.isin(pt, badpt)
        cam, uv = cam[keep], uv[keep]
        _, pt = np.unique(pt[keep], return_inverse=True)
        npt = pt.max() + 1
        X = triangulate(R, t, model, ip, dr, cam, uv, pt, npt)
        say(f"dropped {len(badpt)} tie points with depth outside 0.5-2 x median ({zm:.0f} m)")

    order = np.argsort(pt, kind="stable")          # observations grouped by point (Schur kernel)
    cam, uv, pt = cam[order], uv[order], pt[order]
    ptr = np.r_[0, np.cumsum(np.bincount(pt, minlength=npt))].astype(np.int64)
    schur = None if args.sparse_schur else load_schur()

    # gauge / sanity: GCP reprojection with the ODM solution
    g_res = model.project(cam_xc(R, t, g_cam, g_X), ip, dr) - g_uv
    if len(g_cam):
        say(f"input solution: GCP reprojection RMS {np.sqrt((g_res ** 2).sum(1).mean()):.2f} px")

    ndr = 0 if args.no_dr else model.nq
    NG = 6 * nc + 8 + ndr          # global unknowns: poses, intrinsics, corrections
    o_int, o_dr = 6 * nc, 6 * nc + 8

    def residuals(rv, t, ip, dr, X):
        R = rot(rv)
        r_t = (model.project(cam_xc(R, t, cam, X[pt]), ip, dr) - uv) / SIG_TIE
        r_g = (model.project(cam_xc(R, t, g_cam, g_X), ip, dr) - g_uv) / SIG_GCP
        return r_t, r_g

    def robust_w(r, c):
        e2 = (r ** 2).sum(1) / c ** 2
        return 1 / (1 + e2), np.log1p(e2).sum() * c ** 2

    # constraints that keep the corrections off the smooth (block-bending) modes
    rn = np.arange(model.nr) * RSTEP / 3000.0
    Br = np.stack([rn ** 0, rn, rn ** 3, rn ** 5, rn ** 7], 1)
    gu, gv = np.meshgrid(np.arange(model.gx) / model.gx * 2 - 1, np.arange(model.gy) / model.gy * 2 - 1)
    mono = np.stack([(gu ** i * gv ** j).ravel() for i in range(4) for j in range(4 - i)], 1)   # degree <= 3
    Bfull = np.zeros((model.nq, Br.shape[1] + 2 * mono.shape[1]))
    Bfull[:model.nr, :Br.shape[1]] = Br
    Bfull[model.nr:model.nr + model.ng, Br.shape[1]:Br.shape[1] + mono.shape[1]] = mono
    Bfull[model.nr + model.ng:, Br.shape[1] + mono.shape[1]:] = mono
    Bq, _ = np.linalg.qr(Bfull)
    # smoothness: second differences of dr, and of each grid along x and y
    D = [sparse.diags([1, -2, 1], [0, 1, 2], shape=(model.nr - 2, model.nr)) / SIG_DR2]
    Dx = sparse.kron(sparse.eye(model.gy), sparse.diags([1, -2, 1], [0, 1, 2], shape=(model.gx - 2, model.gx)))
    Dy = sparse.kron(sparse.diags([1, -2, 1], [0, 1, 2], shape=(model.gy - 2, model.gy)), sparse.eye(model.gx))
    for comp in range(2):
        pad = lambda M: sparse.hstack([sparse.csr_matrix((M.shape[0], model.nr + comp * model.ng)), M,
                                       sparse.csr_matrix((M.shape[0], (1 - comp) * model.ng))])
        D += [pad(Dx) / SIG_G2, pad(Dy) / SIG_G2]
    D2s = sparse.vstack([sparse.hstack([D[0], sparse.csr_matrix((D[0].shape[0], 2 * model.ng))])] + D[1:]).tocsr()
    sig0 = np.r_[np.full(model.nr, SIG_DR), np.full(2 * model.ng, SIG_G)]

    def prior_terms(rv, t, dr):
        rows = [((rv - rv0) / SIG_ROT).ravel(), ((t - t0) / SIG_T).ravel()]
        if ndr:
            rows += [dr / sig0, D2s @ dr, Bq.T @ dr / SIG_ORTH]
        if not args.free_intr:
            rows.append((ip - ip0) / SIG_FIX)
        return np.concatenate(rows)

    def cost(rv, t, ip, dr, X, c):
        r_t, r_g = residuals(rv, t, ip, dr, X)
        _, ct = robust_w(r_t, c)
        return ct + (r_g ** 2).sum() + (prior_terms(rv, t, dr) ** 2).sum()

    def jac(rv, t, ip, dr, X, cam_, uvX, fixedX, sig):
        """d(px)/d(params) for observations of cameras cam_ at points uvX.
        returns base (n,2), Jc (n,2,6), Ji (n,2,8), Jd (n,2,2)+node index, Jx (n,2,3).
        Chain rule through Xc = R X + t: the projection is differentiated numerically in Xc only
        (3 projections), then d/dt = dpx/dXc, d/dX = dpx/dXc R, d/drot = dpx/dXc dXc/drot.
        Intrinsics are held (SIG_FIX) unless --free-intr: their columns are then not computed."""
        R = rot(rv)
        Xc = cam_xc(R, t, cam_, uvX)
        base = model.project(Xc, ip, dr)
        n = len(cam_)
        Jxc = np.empty((n, 2, 3))
        for j in range(3):
            d = np.zeros(3); d[j] = 1e-4
            Jxc[:, :, j] = (model.project(Xc + d, ip, dr) - base) / 1e-4
        Jc = np.empty((n, 2, 6)); Ji = np.zeros((n, 2, 8))
        for j in range(3):
            h = 1e-7; d = np.zeros(3); d[j] = h
            Jc[:, :, j] = np.einsum("nck,nk->nc", Jxc, (cam_xc(rot(rv + d), t, cam_, uvX) - Xc) / h)
        Jc[:, :, 3:] = Jxc
        Jx = np.zeros((n, 2, 3)) if fixedX else np.einsum("nck,nkj->ncj", Jxc, R[cam_])
        if args.free_intr:
            for j in range(8):
                h = 1e-7; d = np.zeros(8); d[j] = h
                Ji[:, :, j] = (model.project(Xc, ip + d, dr) - base) / h
        i0, Jd = model.corr_terms(model.brown(Xc, ip))
        return base, Jc / sig, Ji / sig, Jd / sig, i0, Jx / sig

    def sparse_global(cam_, Jc, Ji, Jd, i0, w):
        """(2n x NG) sparse matrix of sqrt(w)-weighted global Jacobian rows."""
        n = len(cam_)
        if n == 0:
            return sparse.csr_matrix((0, NG))
        sw = np.sqrt(w)[:, None, None]
        Jc, Ji, Jd = Jc * sw, Ji * sw, Jd * sw
        rows = np.repeat(np.arange(2 * n), 6 + 8 + (10 if ndr else 0))
        cols_c = cam_[:, None] * 6 + np.arange(6)
        cols_i = np.broadcast_to(o_int + np.arange(8), (n, 8))
        cols = [cols_c, cols_i]
        vals = [Jc, Ji]
        if ndr:
            cols.append(o_dr + i0); vals.append(Jd)
        cols = np.concatenate(cols, 1)                       # (n, k)
        cols = np.repeat(cols, 2, 0).ravel()                 # rows u, v share cols
        vals = np.concatenate(vals, 2).reshape(2 * n, -1).ravel()
        return sparse.csr_matrix((vals, (rows, cols)), shape=(2 * n, NG))

    lam = 1e-3
    c = 3.0          # Cauchy scale, in sigmas
    for it in range(args.iters):
        tic = time.time()
        r_t, r_g = residuals(rv, t, ip, dr, X)
        w_t, _ = robust_w(r_t, c)
        base, Jc, Ji, Jd, i0, Jx = jac(rv, t, ip, dr, X[pt], cam, X[pt], False, SIG_TIE)
        sw = np.sqrt(w_t)
        rt = r_t * sw[:, None]
        Jxw = Jx * sw[:, None, None]
        # GCP rows (fixed points)
        _, gc, gi, gd, g0, _ = jac(rv, t, ip, dr, g_X, g_cam, g_X, True, SIG_GCP)
        Jgcp = sparse_global(g_cam, gc, gi, gd, g0, np.ones(len(g_cam)))
        if schur is not None:
            Gs = [Jc, Ji] + ([Jd] if ndr else [])
            Gk = np.ascontiguousarray(np.concatenate(Gs, 2) * sw[:, None, None])          # (n, 2, K)
            Ck = [cam[:, None] * 6 + np.arange(6),
                  np.broadcast_to(o_int + np.arange(8) if args.free_intr else np.full(8, -1), (len(cam), 8))]
            if ndr:
                Ck.append(o_dr + i0)
            Ck = np.ascontiguousarray(np.concatenate(Ck, 1), dtype=np.int64)
            U = np.zeros((NG, NG)); bg = np.zeros(NG)
            Vinv = np.empty((npt, 3, 3)); bp = np.empty((npt, 3))
            nth = int(max(1, min(os.cpu_count(), 6e9 // (8 * NG * NG))))
            if schur.schur_points(npt, ptr, Gk.shape[2], Ck, Gk, np.ascontiguousarray(Jxw), np.ascontiguousarray(rt),
                                  NG, U, bg, Vinv, bp, nth) != 0:
                raise MemoryError("schur_points: no memory for the thread copies")
            U += (Jgcp.T @ Jgcp).toarray()
            bg += Jgcp.T @ r_g.ravel()

            def wt_dg(dg):        # W^T dg per point = sum_o P_o^T (G_o dg[cols_o])
                q = np.einsum("nck,nk->nc", Gk, np.where(Ck >= 0, dg[Ck], 0))
                out = np.zeros((npt, 3)); np.add.at(out, pt, np.einsum("nci,nc->ni", Jxw, q))
                return out.ravel()
        else:
            Jg = sparse_global(cam, Jc, Ji, Jd, i0, w_t)
            rows = np.repeat(np.arange(2 * len(cam)), 3)
            cols = (np.repeat(pt, 2)[:, None] * 3 + np.arange(3)).ravel()
            Jp = sparse.csr_matrix((Jxw.reshape(-1), (rows, cols)), shape=(2 * len(cam), 3 * npt))
            U = (Jg.T @ Jg).toarray() + (Jgcp.T @ Jgcp).toarray()
            bg = Jg.T @ rt.ravel() + Jgcp.T @ r_g.ravel()
        # priors (linear)
        pr = prior_terms(rv, t, dr)
        dpr = np.zeros(NG)
        Up = np.zeros(NG)
        Up[:6 * nc] = np.tile(np.r_[np.full(3, SIG_ROT ** -2), np.full(3, SIG_T ** -2)], nc)
        dpr[:6 * nc] = np.c_[(rv - rv0) / SIG_ROT ** 2, (t - t0) / SIG_T ** 2].ravel()
        if not args.free_intr:
            Up[o_int:o_int + 8] = SIG_FIX ** -2
            dpr[o_int:o_int + 8] = (ip - ip0) / SIG_FIX ** 2
        U[np.diag_indices(NG)] += Up
        if ndr:
            Pd = np.diag(sig0 ** -2.0) + (D2s.T @ D2s).toarray() + Bq @ Bq.T / SIG_ORTH ** 2
            U[o_dr:, o_dr:] += Pd
            dpr[o_dr:] = Pd @ dr
        bg = bg + dpr
        if schur is not None:
            Sm, rhs = U, bg
        else:
            # 3x3 blocks of Jp^T Jp and W = Jg^T Jp
            Vb = np.zeros((npt, 3, 3))
            np.add.at(Vb, pt, np.einsum("nki,nkj->nij", Jxw, Jxw))
            Vb += 1e-6 * np.eye(3)
            bp = np.zeros((npt, 3)); np.add.at(bp, pt, np.einsum("nki,nk->ni", Jxw, rt))
            W = (Jg.T @ Jp).tocsr()                     # NG x 3npt
            Vinv = np.linalg.inv(Vb)
            bc = np.broadcast_to(np.arange(npt)[:, None, None] * 3 + np.arange(3)[None, None, :], (npt, 3, 3))
            Vi = sparse.csr_matrix((Vinv.ravel(), (np.repeat(np.arange(3 * npt), 3), bc.ravel())),
                                   shape=(3 * npt, 3 * npt))
            WVi = W @ Vi
            Sm = U - (WVi @ W.T).toarray()
            rhs = bg - WVi @ bp.ravel()

            def wt_dg(dg):
                return W.T @ dg
        ok = False
        c0 = cost(rv, t, ip, dr, X, c)
        for _ in range(8):
            A = Sm + lam * np.diag(np.maximum(np.diag(Sm), 1e-9))
            dg = -np.linalg.solve(A, rhs)
            dp = -(Vinv @ (bp + wt_dg(dg).reshape(npt, 3))[..., None])[..., 0]
            nrv = rv + dg[:6 * nc].reshape(nc, 6)[:, :3]
            nt = t + dg[:6 * nc].reshape(nc, 6)[:, 3:]
            nip = ip + dg[o_int:o_int + 8]
            ndr_ = dr + (dg[o_dr:] if ndr else 0)
            nX = X + dp
            c1 = cost(nrv, nt, nip, ndr_, nX, c)
            if c1 < c0:
                rv, t, ip, dr, X = nrv, nt, nip, ndr_, nX
                lam = max(lam / 3, 1e-7); ok = True
                break
            lam *= 10
        r_t, r_g = residuals(rv, t, ip, dr, X)
        e = np.hypot(*r_t.T) * SIG_TIE
        say(f"it {it}: cost {c0:.4g} -> {c1:.4g} lam {lam:.1e} tie |res| median {np.median(e):.3f} px "
            f"p90 {np.percentile(e, 90):.3f} GCP RMS {np.sqrt((r_g ** 2).sum(1).mean()) * SIG_GCP:.2f} px "
            f"f {ip[0]:.6f} dr max {np.abs(dr).max():.3f} px  {time.time() - tic:.0f}s")
        if not ok or abs(c0 - c1) / c0 < args.tol:
            break

    # ---------------- report
    R = rot(rv)
    C0 = -np.einsum("nji,nj->ni", rot(rv0), t0); C1 = -np.einsum("nji,nj->ni", R, t)
    say("camera centre shift vs ODM (m):")
    for k, n in enumerate(names):
        dC = C1[k] - C0[k]
        say(f"  {n} dx {dC[0]:+.3f} dy {dC[1]:+.3f} dz {dC[2]:+.3f} rot {np.degrees(np.linalg.norm(rv[k]-rv0[k]))*3600:.0f}\"")
    say("intrinsics ODM -> BBA: " + ", ".join(f"{k} {a:.6g}->{b:.6g}" for k, a, b in zip(INTR, ip0, ip)))
    say("dr nodes (px): " + " ".join(f"{x:+.3f}" for x in dr[:model.nr]))
    gq = dr[model.nr:].reshape(2, model.gy, model.gx)
    say(f"grid correction: x std {gq[0].std():.3f} max {np.abs(gq[0]).max():.3f}, y std {gq[1].std():.3f} max {np.abs(gq[1]).max():.3f} px")
    gres = model.project(cam_xc(R, t, g_cam, g_X), ip, dr) - g_uv if len(g_cam) else np.zeros((0, 2))
    gid = np.array(g_id); ug = sorted(set(g_id)) if len(g_cam) else []
    gi = np.array([ug.index(x) for x in g_id])
    Xg = triangulate(R, t, model, ip, dr, g_cam, g_uv, gi, len(ug)) if len(ug) else []
    for k, u in enumerate(ug):
        Xk = g_X[gid == u][0]
        say(f"  GCP {u} triangulated - surveyed: dX {Xg[k, 0] - Xk[0]:+.3f} dY {Xg[k, 1] - Xk[1]:+.3f} dZ {Xg[k, 2] - Xk[2]:+.3f} m")
    for k in range(len(g_cam)):
        say(f"  GCP {g_id[k]} {names[g_cam[k]]} res {gres[k, 0]:+.2f} {gres[k, 1]:+.2f} px")

    # ---------------- write
    newcam = dict(camd)
    for k, v in zip(INTR, ip):
        newcam[k] = float(v)
    newcam["focal_y"] = float(ip[0])
    newcam["kjr_dr_step"] = RSTEP
    newcam["kjr_dr"] = [float(x) for x in dr[:model.nr]]
    newcam["kjr_grid_step"] = GSTEP
    newcam["kjr_grid_shape"] = [model.gy, model.gx]
    newcam["kjr_grid_x"] = [float(x) for x in dr[model.nr:model.nr + model.ng]]
    newcam["kjr_grid_y"] = [float(x) for x in dr[model.nr + model.ng:]]
    shots = {}
    for k, n in enumerate(names):
        s = dict(rec["shots"][n]); s["rotation"] = [float(x) for x in rv[k]]; s["translation"] = [float(x) for x in t[k]]
        shots[n] = s
    json.dump({"intrinsics": newcam, "shots": {n: {"rotation": s["rotation"], "translation": s["translation"]} for n, s in shots.items()}},
              open(os.path.join(args.out, "bba.json"), "w"), indent=1)
    od = os.path.join(args.out, "odm")
    os.makedirs(os.path.join(od, "opensfm/undistorted"), exist_ok=True)
    os.makedirs(os.path.join(od, "odm_georeferencing"), exist_ok=True)
    shutil.copy(os.path.join(args.odm, "odm_georeferencing/coords.txt"), os.path.join(od, "odm_georeferencing/coords.txt"))
    if not os.path.exists(os.path.join(od, "odm_orthophoto")):      # output grid reference only (read only)
        os.symlink(os.path.join(os.path.abspath(args.odm), "odm_orthophoto"), os.path.join(od, "odm_orthophoto"))
    np.save(os.path.join(args.out, "points.npy"), X.astype(np.float32))   # all tie points (kjr/refdsm.py)
    pts = {}   # a sparse subset of the adjusted tie points (ground plane fit in stereo.py)
    sel = np.linspace(0, npt - 1, min(npt, 100000)).astype(int)
    for k in sel:
        pts[str(k)] = {"coordinates": [float(x) for x in X[k]], "color": [0, 0, 0]}
    json.dump([{"cameras": {cam_key: newcam}, "shots": shots, "points": pts,
                "reference_lla": rec.get("reference_lla", {})}],
              open(os.path.join(od, "opensfm/reconstruction.json"), "w"))
    ucam = {"projection_type": "perspective", "width": camd["width"], "height": camd["height"],
            "focal": float(ip[0]), "k1": 0.0, "k2": 0.0}
    ushots = {}
    for n, s in shots.items():
        u = dict(s); ushots[n + ".tif"] = u
    json.dump([{"cameras": {cam_key: ucam}, "shots": ushots, "points": pts}],
              open(os.path.join(od, "opensfm/undistorted/reconstruction.json"), "w"))
    say(f"written {args.out}/bba.json and {od}")


if __name__ == "__main__":
    main()
