#!/usr/bin/env bash
# uav_kjr main entry: one block from the raw images and a GCP list, no ODM.
#   1. upstream OpenSfM (BSD) SfM with the GCPs     tools/run_opensfm.sh -> <work>/sfm/{osfm,odm}
#   2. ties, own BBA, stereo, fusion, blend, ortho, products   run_block.sh
# Usage: run.sh <BLOCK> <image dir> <gcp list> <work dir> <products dir>
#   gcp list: OpenSfM/ODM format (first line projection, then E N Z px py image [id])
# Env: OSFM_ENV, OSFM_SRC (tools/build_opensfm.sh); GSD=<m> ortho cell (default flying height / focal);
#      STAGES, SIGMA, BBA_ARGS ... as run_block.sh
set -euo pipefail
[ $# -eq 5 ] || { sed -n 2,8p "$0"; exit 1; }
B="$1"; SRC="$2"; GCP="$3"; W=$(realpath -m "$4"); OUT="$5"
HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$W"
[ -f "$W/sfm/odm/odm_dem/dsm.tif" ] || GSD=${GSD:-} "$HERE/tools/run_opensfm.sh" "$SRC" "$GCP" "$W/sfm" 2>&1 | tee -a "$W/sfm.log"
REF_DSM="$W/sfm/odm/odm_dem/dsm.tif" "$HERE/run_block.sh" "$B" "$W/sfm/odm" "$SRC" "$GCP" "$W" "$OUT"
