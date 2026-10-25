#!/usr/bin/env python3
"""ODM-KJR BBA step 1: sub-pixel tie points on the ORIGINAL (distorted) images.

UAV library tiepoints.match_pair: SIFT -> MAGSAC++ (F, 1 px) -> Gruen ALSC (cuda).
Pairs: the overlapping pairs of a stereo run (stereo_meta.json), which include
the cross-strip pairs. Each pair keeps its own 2-view ties (ALSC moves only the
second image's point; multi-view tracks are grouped on the first image's unmoved
SIFT position, kjr/bba.py).
Writes <out>/ties/<a>__<b>.npz with pa, pb (n,2) px (x = column, y = row).
Usage: ties.py <stereo meta json> <original image dir> <out> [--workers 3]
"""
import argparse, json, os, sys, time
import multiprocessing as mp
import cv2
import numpy as np

MAXP = None     # keep all: tracks are grouped on the first image's exact SIFT
                # position, so every pair must keep the same anchors (no subsampling)


def grey(path):
    im = cv2.imread(path, cv2.IMREAD_COLOR).astype(np.float32)
    return im @ np.array([0.114, 0.587, 0.299], np.float32)


def job(t):
    a, b, src, out = t
    dst = os.path.join(out, f"{a}__{b}.npz")
    if os.path.exists(dst):
        return a, b, None
    from tiepoints.tiepoints import match_pair
    t0 = time.time()
    A, B = grey(os.path.join(src, a)), grey(os.path.join(src, b))
    for k in range(20):
        try:
            tp = match_pair(A, B, backend="sift", alsc="cuda", max_points=MAXP)
            break
        except RuntimeError:          # ALSC status -100: cudaMalloc failed, another
            time.sleep(5 + 10 * np.random.rand())   # worker holds the 8 GB GPU; retry
    else:
        raise RuntimeError(f"{a} {b}: ALSC kept failing")
    np.savez(dst, pa=tp.pi.astype(np.float64), pb=tp.pj.astype(np.float64))
    return a, b, (tp.stats.get("ransac_F", 0), tp.stats.get("final", 0), time.time() - t0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("meta"); ap.add_argument("src"); ap.add_argument("out")
    ap.add_argument("--workers", type=int, default=3)
    a = ap.parse_args()
    os.makedirs(os.path.join(a.out, "ties"), exist_ok=True)
    pairs = [(p[0].replace(".tif", ""), p[1].replace(".tif", "")) for p in json.load(open(a.meta))["pairs"]]
    jobs = [(x, y, a.src, os.path.join(a.out, "ties")) for x, y in pairs]
    with mp.get_context("spawn").Pool(a.workers) as pool:   # fork breaks the CUDA ALSC
        for x, y, s in pool.imap_unordered(job, jobs):
            print(x, y, "skip" if s is None else f"F {s[0]} kept {s[1]} {s[2]:.0f}s", flush=True)


if __name__ == "__main__":
    main()
