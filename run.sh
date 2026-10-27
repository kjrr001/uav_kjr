#!/usr/bin/env bash
# Old command line, kept for compatibility: writes a control file and runs kjr_run (use kjr_run).
# Usage: run.sh <BLOCK> <image dir> <gcp list> <work dir> <products dir>
# Env: GSD, SIGMA, BBA_ITERS, TIE_WORKERS, BBA_ARGS, SHAPE_DSM, REF_DSM; EXISTING (default skip), DRY_RUN=1;
#      START / END (stage names, kjr_run --stages)
set -euo pipefail
[ $# -eq 5 ] || { sed -n 2,6p "$0"; exit 1; }
HERE=$(cd "$(dirname "$0")" && pwd); W=$(realpath -m "$4"); mkdir -p "$W/job"
C="$W/job/run_sh.cfg"
cat > "$C" <<CFG
BLOCK        : $1
INPUT_DIR    : $(realpath "$2")
GCP_FILE     : $(realpath "$3")
WORK_DIR     : $W
OUTPUT_DIR   : $(realpath -m "$5")
SFM_SOURCE   : opensfm
EXISTING     : ${EXISTING:-skip}
START        : ${START:-}
END          : ${END:-DIAG}
GSD          : ${GSD:-auto}
REF_DSM      : ${REF_DSM:-auto}
SHAPE_DSM    : ${SHAPE_DSM:-ref}
BLEND_SIGMA  : ${SIGMA:-30}
BBA_ITERS    : ${BBA_ITERS:-30}
TIE_WORKERS  : ${TIE_WORKERS:-2}
BBA_OPTIONS  : ${BBA_ARGS:-}
CFG
exec "$HERE/kjr_run" "$C" ${DRY_RUN:+--dry-run}
