#!/usr/bin/env bash
# ODM-KJR: one block, end to end, from a finished ODM GCP run (read only).
#
#   ref      stereo + fuse with the ODM cameras    -> work/ref   (pair list, long-wave reference DSM)
#   ties     sub-pixel ties on the original JPEGs  -> work/ties  (SIFT -> MAGSAC++ -> ALSC, GPU)
#   bba      own bundle adjustment + radial/grid image correction -> work/bba (+ work/bba/odm)
#   undist   Lanczos undistortion with the new camera            -> work/images
#   stereo   census-SGM per pair on the new images, search centred on the ODM DSM -> work/stereo/pairs
#   fuse     robust fusion of the pair rasters                    -> work/stereo/dsm_fused.tif
#   blend    long waves (> sigma) from ref, detail from stereo    -> work/final/dsm_fused.tif
#   ortho    orthophoto from the final DSM                        -> work/final/ortho.tif
#   products COG / ENVI / JPG deliverables                        -> <products dir>
#   diag     colour hillshade, GCP heights, ring test             -> work/diag
#
# Usage: run_block.sh <BLOCK> <odm run> <image dir> <gcp list> <work dir> <products dir>
# Env:   STAGES="ties bba ..." run only these (default: all, in order)
#        REF_DSM=<tif>  use this as the long-wave reference instead of making ref/dsm_fused.tif
#                       (recommended: <odm run>/odm_dem/dsm.tif; JJ-1 GCP RMS 8.3 vs 13.6 cm, 45 min less)
#        SIGMA=30 (blend, m)  BBA_ITERS=30  TIE_WORKERS=2
# Every stage skips when its output already exists: rerun after a failure resumes;
# delete a stage's output to redo it.
set -euo pipefail
[ $# -eq 6 ] || { sed -n 2,25p "$0"; exit 1; }
B="$1"; ODM=$(realpath "$2"); SRC=$(realpath "$3"); GCP=$(realpath "$4"); W="$5"; OUT="$6"
HERE=$(cd "$(dirname "$0")" && pwd)
source "$HERE/env.sh"
SIGMA=${SIGMA:-30}; BBA_ITERS=${BBA_ITERS:-30}; TIE_WORKERS=${TIE_WORKERS:-2}
STAGES=${STAGES:-"ref ties bba undist stereo fuse blend ortho products diag"}
mkdir -p "$W" "$OUT"; W=$(realpath "$W")
LOG="$W/run_block.log"
F='Warn|warn|census variant|exact streaming|FutureWarning|diags_array|warnings.warn|^Note: In Python'
K="$HERE/kjr"; T="$HERE/tools"
want() { [[ " $STAGES " == *" $1 "* ]]; }
say() { echo "== $(date '+%F %T') $B $*" | tee -a "$LOG"; }
run() {   # run a step, log its output without the noise, stop on failure
  set +e; "$@" 2>&1 | grep -vE "$F" >> "$LOG"; local s=${PIPESTATUS[0]}; set -e
  [ "$s" -eq 0 ] || { say "FAILED ($s): $*"; exit 1; }
}

for f in "$ODM/opensfm/reconstruction.json" "$ODM/opensfm/undistorted/reconstruction.json" \
         "$ODM/odm_georeferencing/coords.txt" "$ODM/odm_orthophoto/odm_orthophoto.tif" "$GCP"; do
  [ -e "$f" ] || { echo "missing: $f"; exit 1; }
done
say "start: odm $ODM, images $SRC, gcp $GCP, stages: $STAGES"
REF=${REF_DSM:-$W/ref/dsm_fused.tif}
PRIOR=${REF_DSM:-$ODM/odm_dem/dsm.tif}     # centres the stereo search (+-8 px = only ~+-1 m of height)

if want ref; then
  if [ -n "${REF_DSM:-}" ]; then
    [ -f "$W/ref/stereo_meta.json" ] || { say "ref: pair list only (REF_DSM given)"; run $PY "$K/stereo.py" "$ODM" "$W/ref" --meta-only --prior-dsm "$PRIOR"; }
  elif [ ! -f "$W/ref/dsm_fused.tif" ]; then
    say "ref stereo (ODM cameras)"; run $PY "$K/stereo.py" "$ODM" "$W/ref" --prior-dsm "$PRIOR"
    say "ref fuse";                 run $PY "$K/fuse.py" "$W/ref"
  fi
fi
if want ties; then
  N=$($PY -c "import json; print(len(json.load(open('$W/ref/stereo_meta.json'))['pairs']))")
  if [ "$(ls "$W/ties" 2>/dev/null | grep -c npz)" -lt "$N" ]; then
    say "ties ($N pairs)"; run $PY "$K/ties.py" "$W/ref/stereo_meta.json" "$SRC" "$W" --workers "$TIE_WORKERS"
  fi
fi
if want bba && [ ! -f "$W/bba/bba.json" ]; then
  say "bba ($BBA_ITERS iterations)"; run $PY "$K/bba.py" "$ODM" "$W/ties" "$GCP" "$W/bba" --iters "$BBA_ITERS"
  grep -E "^it |GCP . triangulated|grid correction" "$W/bba/bba.log" | tail -8 | tee -a "$LOG"
fi
if want undist; then
  NAMES=$(ls "$SRC" | grep -i '\.jpg$' | sed 's/$/.tif/')
  say "undistort"; run $PY "$K/undistort.py" "$W/bba/odm" "$SRC" "$W/images" $NAMES
fi
if want stereo; then say "stereo"; run $PY "$K/stereo.py" "$W/bba/odm" "$W/stereo" --img-dir "$W/images" --prior-dsm "$PRIOR"; fi
if want fuse && [ ! -f "$W/stereo/dsm_fused.tif" ]; then say "fuse"; run $PY "$K/fuse.py" "$W/stereo"; fi
if want blend && [ ! -f "$W/final/dsm_fused.tif" ]; then
  mkdir -p "$W/final"
  say "blend (sigma $SIGMA m, reference $REF)"; run $PY "$K/lowblend.py" "$W/stereo/dsm_fused.tif" "$REF" "$W/final/dsm_fused.tif" --sigma "$SIGMA"
fi
if want ortho && [ ! -f "$W/final/ortho.tif" ]; then
  say "ortho"; run $PY "$K/ortho.py" "$W/bba/odm" "$W/final" --img-dir "$W/images"
fi
if want products; then say "products -> $OUT"; run "$K/make_products.sh" "$B" "$W/final" "$OUT"; fi
if want diag; then
  mkdir -p "$W/diag"; say "diag"
  run $PY "$T/color_hillshade.py" "$W/final/dsm_fused.tif" "$OUT/${B}_kjr_dsm_color_hillshade_small.jpg"
  $PY "$T/gcp_check.py" "$GCP" "$REF" "$W/stereo/dsm_fused.tif" "$W/final/dsm_fused.tif" 2>&1 | grep -vE "$F" | tee -a "$LOG"
  $PY "$T/diag_imgfixed.py" "$W/final/dsm_fused.tif" "$W/bba/odm" "$W/diag/imgfixed.npy" 2>&1 | grep "hp DSM" | tee -a "$LOG"
fi
say "done"
