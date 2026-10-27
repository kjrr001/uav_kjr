#!/usr/bin/env bash
# Deliverable folder: per block DSM, ortho, small colour hillshade, small ortho (same size), the
# camera-image footprint map, plus one context map of all blocks on a regional DEM.
# (one block at a time: kjr_run with EXPORT_DIR, END EXPORT)
# Usage: export_products.sh <products root> <work root> <out root> <regional dem.tif> BLOCK [BLOCK ...]
#   <products root>/<B>/<B>_kjr_{dsm,ortho}.tif, _dsm_color_hillshade_small.jpg; <work root>/<B>/bba/odm
set -euo pipefail
HERE=$(cd "$(dirname "$0")/.." && pwd); source "$HERE/env.sh"
SRC="$1"; WRK="$2"; DST="$3"; DEM="$4"; shift 4
CTX=()
for B in "$@"; do
  mkdir -p "$DST/$B"
  $PY "$HERE/tools/small_ortho.py" "$SRC/$B/${B}_kjr_ortho.tif" "$SRC/$B/${B}_kjr_dsm_color_hillshade_small.jpg" "$DST/$B/${B}_kjr_ortho_small.jpg"
  cp -v "$SRC/$B/${B}_kjr_dsm.tif" "$SRC/$B/${B}_kjr_ortho.tif" "$SRC/$B/${B}_kjr_dsm_color_hillshade_small.jpg" "$DST/$B/"
  $PY "$HERE/tools/image_footprints.py" "${B}_kjr" "$WRK/$B" "$SRC/$B" "$DST/$B/${B}_image_footprints.png" "$DST/$B/${B}_kjr_ortho_small.jpg"
  CTX+=("$B=$SRC/$B/${B}_kjr_dsm.tif")
done
$PY "$HERE/tools/context_map.py" "$DEM" "$DST/uav_context_map.png" "${CTX[@]}"
