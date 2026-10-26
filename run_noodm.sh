#!/usr/bin/env bash
# ODM-KJR without ODM: one block from the raw DJI JPEGs and a GCP list.
#
#   init     cameras from EXIF/XMP (GPS, gimbal yaw/pitch, DJI calibrated focal) -> work/init
#   pairs    overlapping pairs from the EXIF footprints                          -> work/ref/stereo_meta.json
#   ties, tracks                      (as run_block.sh)                          -> work/ties, work/tracks.npz
#   rel      bundle adjustment on the ties, GPS as a weak prior (5 m, 0.1 rad), no GCP
#                                                                                -> work/rel/odm
#   abs      GCPs triangulated from it -> similarity onto the survey (datum, heights), then
#            bundle adjustment with the GCPs and self-calibration (f, c, k1-k3, p1, p2), the GCPs
#            weighted as sigma SIG_GCP (default 0.03 px): at the marking sigma (0.3 px) 20 GCP
#            observations are outweighed by millions of ties and the radial terms bend the block
#            (KJ-3: GCP heights -0.3..+0.7 m, tie-point DSM vs ODM NMAD 57 cm -> 13 cm at 0.03 px)
#                                                                                -> work/abs/odm
#   refdsm   reference DSM from the adjusted tie points (2 m) -> work/abs/odm/odm_dem/dsm.tif
#            (centres the stereo search only: as long-wave shape it made KJ-3 worse, GCP RMS
#            6.1 -> 13.2 cm, so blend keeps the stereo DSM's own shape: SHAPE_DSM=self)
#   then run_block.sh on work/abs/odm (bba with the image corrections, stereo, fuse, blend,
#   ortho, products, diag), REF_DSM = that reference DSM.
#
# Usage: run_noodm.sh <BLOCK> <image dir> <gcp list> <work dir> <products dir>
# Env:   CAMERA=<calibration json>  OpenSfM/ODM cameras.json / reconstruction.json or kjr bba.json of the
#                 same camera (calibration flight, earlier project): intrinsics held at it. Without
#                 it they are self-calibrated, which nadir parallel strips + few GCPs determine
#                 badly (KJ-3: check points RMS 80 cm against 21 cm with a calibration)
#        CHECK=1  leave-one-out GCP check points after the adjustment (+5-10 min per GCP)
#        SIG_GCP=0.03 (px)  EXIF_ARGS="--agl 120" (no RelativeAltitude tag)  STAGES (run_block part)
# Needs: one camera model, near-nadir images with GPS, >= 3 GCPs seen in >= 2 images each.
set -euo pipefail
[ $# -eq 5 ] || { sed -n 2,21p "$0"; exit 1; }
B="$1"; SRC=$(realpath "$2"); GCP=$(realpath "$3"); W="$4"; OUT="$5"
HERE=$(cd "$(dirname "$0")" && pwd)
source "$HERE/env.sh"
SIG_GCP=${SIG_GCP:-0.03}
K="$HERE/kjr"; mkdir -p "$W"; W=$(realpath "$W"); LOG="$W/run_noodm.log"
say() { echo "== $(date '+%F %T') $B $*" | tee -a "$LOG"; }
run() { set +e; "$@" 2>&1 | grep -vE 'Warn|warn|census variant' >> "$LOG"; local s=${PIPESTATUS[0]}; set -e
        [ "$s" -eq 0 ] || { say "FAILED ($s): $*"; exit 1; }; }

CAM_ARGS=(); [ -n "${CAMERA:-}" ] && CAM_ARGS=(--camera "$(realpath "$CAMERA")")
[ -f "$W/init/opensfm/reconstruction.json" ] || { say "init (EXIF${CAMERA:+, camera $CAMERA})"; run $PY "$K/exif_init.py" "$SRC" "$W/init" "${CAM_ARGS[@]}" ${EXIF_ARGS:-}; }
# intrinsics: held at the given calibration, else self-calibrated (weakly determined, see MANUAL 4.1)
if [ -n "${CAMERA:-}" ]; then INTR=(); else INTR=(--free-intr); say "no CAMERA given: self-calibration - check the result with CHECK=1"; fi
[ -f "$W/ref/stereo_meta.json" ] || { say "pairs"; run $PY "$K/stereo.py" "$W/init" "$W/ref" --meta-only --prior-dsm none --min-overlap 0.05; }
N=$($PY -c "import json; print(len(json.load(open('$W/ref/stereo_meta.json'))['pairs']))")
if [ "$(ls "$W/ties" 2>/dev/null | grep -c npz)" -lt "$N" ]; then
  say "ties ($N pairs)"; run $PY "$K/ties.py" "$W/ref/stereo_meta.json" "$SRC" "$W" --workers 2 --keep-sift
fi
[ -f "$W/tracks.npz" ] || { say "tracks"; run $PY "$K/tracks.py" "$W/ref/stereo_meta.json" "$SRC" "$W"; }
rm -rf "$W/sift"
[ -f "$W/rel/bba.json" ] || { say "rel (ties + GPS)"; run $PY "$K/bba.py" "$W/init" "$W/ties" "$GCP" "$W/rel" \
    --tracks "$W/tracks.npz" --no-gcp --no-dr --sig-t 5 --sig-rot 0.1 --iters 30; }
[ -f "$W/abs/bba.json" ] || { say "abs (GCP alignment + self-calibration)"; run $PY "$K/bba.py" "$W/rel/odm" "$W/ties" "$GCP" "$W/abs" \
    --tracks "$W/tracks.npz" --align-gcp "${INTR[@]}" --no-dr --sig-t 50 --sig-rot 0.5 --sig-gcp "$SIG_GCP" --iters 30; }
grep -E "align-gcp|^it |GCP . triangulated|intrinsics" "$W/abs/bba.log" | tail -12 | tee -a "$LOG"
[ -f "$W/abs/odm/odm_dem/dsm.tif" ] || { say "refdsm"; run $PY "$K/refdsm.py" "$W/abs" "$W/init" --out "$W/abs/odm/odm_dem/dsm.tif"; }
cp -n "$W/init/kjr_exif.json" "$W/abs/odm/" 2>/dev/null || true
if [ "${CHECK:-0}" = 1 ] && [ ! -f "$W/check_points.txt" ]; then      # leave-one-out GCP check points
  say "check points (each GCP left out of the abs adjustment in turn)"
  for h in $(awk 'NR>1 && NF>=6 {print (NF>6 ? $7 : $6)}' "$GCP" | sort -u); do
    run $PY "$K/bba.py" "$W/rel/odm" "$W/ties" "$GCP" "$W/check/$h" --tracks "$W/tracks.npz" --align-gcp "${INTR[@]}" \
        --no-dr --sig-t 50 --sig-rot 0.5 --sig-gcp "$SIG_GCP" --iters 30 --holdout "$h"
    grep "CHECK POINT" "$W/check/$h/bba.log" >> "$W/check_points.txt" || true
  done
  $PY -c "
import re, numpy as np
d = [float(re.search(r'dZ ([-+.0-9]+)', l).group(1)) for l in open('$W/check_points.txt')]
print(f'check points: {len(d)}, dZ RMS {np.sqrt(np.mean(np.square(d)))*100:.1f} cm, max |dZ| {np.max(np.abs(d))*100:.1f} cm')" | tee -a "$LOG"
  rm -rf "$W/check"
fi
say "run_block on $W/abs/odm"
BBA_ARGS="--sig-gcp $SIG_GCP" SHAPE_DSM=self REF_DSM="$W/abs/odm/odm_dem/dsm.tif" STAGES="${STAGES:-bba undist stereo fuse blend ortho products diag}" \
  "$HERE/run_block.sh" "$B" "$W/abs/odm" "$SRC" "$GCP" "$W" "$OUT"
say "done"
