#!/usr/bin/env bash
# Deliverable folder: per block DSM, ortho, small colour hillshade, small ortho (same size), the
# camera-image footprint map, plus one context map of all blocks on a regional DEM.
# Usage: export_products.sh <products root> <work root> <out root> <regional dem.tif> BLOCK [BLOCK ...]
#   <products root>/<B>/<B>_kjr_{dsm,ortho}.tif, _dsm_color_hillshade_small.jpg; <work root>/<B>/bba/odm
set -euo pipefail
HERE=$(cd "$(dirname "$0")/.." && pwd); source "$HERE/env.sh"
SRC="$1"; WRK="$2"; DST="$3"; DEM="$4"; shift 4
CTX=()
for B in "$@"; do
  mkdir -p "$DST/$B"
  $PY - "$SRC/$B" "$DST/$B" "$B" <<'PY'
import sys
import numpy as np
from osgeo import gdal
from PIL import Image
gdal.UseExceptions()
s, d, b = sys.argv[1:]
w, h = Image.open(f"{s}/{b}_kjr_dsm_color_hillshade_small.jpg").size
a = gdal.Translate("", f"{s}/{b}_kjr_ortho.tif", format="MEM", width=w, height=h, resampleAlg="average").ReadAsArray()
rgb = np.where(a[3:4] > 127, a[:3], 255).transpose(1, 2, 0).astype(np.uint8)
Image.fromarray(rgb).save(f"{d}/{b}_kjr_ortho_small.jpg", quality=90)
PY
  cp -v "$SRC/$B/${B}_kjr_dsm.tif" "$SRC/$B/${B}_kjr_ortho.tif" "$SRC/$B/${B}_kjr_dsm_color_hillshade_small.jpg" "$DST/$B/"
  $PY "$HERE/tools/image_footprints.py" "$B" "$WRK/$B" "$SRC/$B" "$DST/$B/${B}_image_footprints.png" "$DST/$B/${B}_kjr_ortho_small.jpg"
  CTX+=("$B=$SRC/$B/${B}_kjr_dsm.tif")
done
$PY "$HERE/tools/context_map.py" "$DEM" "$DST/uav_context_map.png" "${CTX[@]}"
