#!/usr/bin/env bash
# Old command line, kept for compatibility: one block from a finished ODM GCP run (read only).
# Writes a control file with SFM_SOURCE = the ODM run and runs kjr_run (use kjr_run).
# Usage: run_block.sh <BLOCK> <odm run> <image dir> <gcp list> <work dir> <products dir>
# Env:   STAGES="ties bba ..." -> START/END = first/last of them (ref ties tracks bba undist stereo
#        fuse blend ortho products diag; naming tracks turns USE_TRACKS on)
#        REF_DSM=<tif> (unset: reference DSM by stereo with the ODM cameras, as before), SHAPE_DSM,
#        SIGMA, BBA_ITERS, TIE_WORKERS, BBA_ARGS; EXISTING (default skip), DRY_RUN=1
set -euo pipefail
[ $# -eq 6 ] || { sed -n 2,9p "$0"; exit 1; }
HERE=$(cd "$(dirname "$0")" && pwd); W=$(realpath -m "$5"); mkdir -p "$W/job"
declare -A M=([ref]=PAIRS [ties]=TIEPOINT [tracks]=TRACKS [bba]=BBA [undist]=UNDISTORT [stereo]=DSM_EXTRACTION
              [fuse]=FUSE [blend]=BLEND [ortho]=ORTHO [products]=PRODUCTS [diag]=DIAG)
ST=(${STAGES:-ref diag}); TRK=no
for s in "${ST[@]}"; do [ -n "${M[$s]:-}" ] || { echo "unknown stage $s"; exit 1; }; [ "$s" = tracks ] && TRK=yes; done
C="$W/job/run_block.cfg"
cat > "$C" <<CFG
BLOCK        : $1
SFM_SOURCE   : $(realpath "$2")
INPUT_DIR    : $(realpath "$3")
GCP_FILE     : $(realpath "$4")
WORK_DIR     : $W
OUTPUT_DIR   : $(realpath -m "$6")
EXISTING     : ${EXISTING:-skip}
START        : ${M[${ST[0]}]}
END          : ${M[${ST[-1]}]}
USE_TRACKS   : $TRK
REF_DSM      : ${REF_DSM:-stereo}
SHAPE_DSM    : ${SHAPE_DSM:-ref}
BLEND_SIGMA  : ${SIGMA:-30}
BBA_ITERS    : ${BBA_ITERS:-30}
TIE_WORKERS  : ${TIE_WORKERS:-2}
BBA_OPTIONS  : ${BBA_ARGS:-}
CFG
exec "$HERE/kjr_run" "$C" ${DRY_RUN:+--dry-run}
