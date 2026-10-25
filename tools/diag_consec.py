"""Pose consistency along a strip: pairs (i-1,i) and (i,i+1) share image i;
median height difference where their rasters overlap (after the fuse ERODE).
A pose error shows as a step of several cm; KJ-1 v3 had up to 18 cm, the final run < 8 cm.
Image names DJI_%04d.JPG (DJI default). Usage: diag_consec.py <stereo dir> <first> <last>   e.g. 53 64"""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from kjr.fuse import load_pair
d, a, b = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
nm = lambda k: f"DJI_{k:04d}.JPG.tif"
out = []
for i in range(a + 1, b):
    f1 = os.path.join(d, "pairs", f"{nm(i-1)}__{nm(i)}.npz"); f2 = os.path.join(d, "pairs", f"{nm(i)}__{nm(i+1)}.npz")
    if not (os.path.exists(f1) and os.path.exists(f2)):
        out.append("  n/a"); continue
    r1, c1, z1, _ = load_pair(f1, None); r2, c2, z2, _ = load_pair(f2, None)
    R0, C0 = max(r1, r2), max(c1, c2); R1 = min(r1 + z1.shape[0], r2 + z2.shape[0]); C1 = min(c1 + z1.shape[1], c2 + z2.shape[1])
    if R1 <= R0 or C1 <= C0:
        out.append("  n/a"); continue
    dz = z2[R0 - r2:R1 - r2, C0 - c2:C1 - c2] - z1[R0 - r1:R1 - r1, C0 - c1:C1 - c1]
    dz = dz[np.isfinite(dz)]
    out.append(f"{np.median(dz) * 100:+5.1f}" if len(dz) > 1000 else "  n/a")
print(f"strip {a}-{b}:", " ".join(out), "cm")
