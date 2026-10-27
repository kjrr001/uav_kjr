#!/usr/bin/env bash
# Upstream OpenSfM (mapillary/OpenSfM, BSD-2) on one block, then the kjr run dir from it
# (kjr/osfm_import.py). Replaces ODM's SfM for run_block.sh.
#   camera model: brown (k1-k3, p1, p2) started from the EXIF focal, self-calibrated with the GCPs
# Usage: run_opensfm.sh <image dir> <gcp list (ODM/OpenSfM format)> <work dir>
#   -> <work dir>/osfm (OpenSfM dataset), <work dir>/odm (run dir for run_block.sh)
# Env: OSFM_ENV  conda prefix with OpenSfM built (default /work/eridu_drone/kjr/opensfm/env)
#      OSFM_SRC  OpenSfM source checkout (bin/opensfm_main.py)
set -euo pipefail
[ $# -eq 3 ] || { sed -n 2,9p "$0"; exit 1; }
SRC=$(realpath "$1"); GCP=$(realpath "$2"); W=$(realpath -m "$3")
HERE=$(cd "$(dirname "$0")/.." && pwd)
OSFM_ENV=${OSFM_ENV:-/work/eridu_drone/kjr/opensfm/env}; OSFM_SRC=${OSFM_SRC:-/work/eridu_drone/kjr/opensfm/src}
D=$W/osfm; mkdir -p "$D/images"
for f in "$SRC"/*; do case "${f,,}" in *.jpg|*.jpeg|*.tif|*.tiff|*.png) ln -sf "$f" "$D/images/";; esac; done
cp "$GCP" "$D/gcp_list.txt"
{ cat "$HERE/tools/opensfm_config.yaml"; echo "processes: $(( $(nproc) > 2 ? $(nproc) - 2 : 1 ))"; } > "$D/config.yaml"
osfm() { echo "== $(date '+%F %T') opensfm $1"; PATH=$OSFM_ENV/bin:$PATH python "$OSFM_SRC/bin/opensfm_main.py" "$1" "$D" > "$D/$1.log" 2>&1 \
         || { echo "opensfm $1 FAILED:"; tail -20 "$D/$1.log"; exit 1; }; }
if [ ! -f "$D/camera_models_overrides.json" ]; then        # every EXIF camera as brown, same focal
  osfm extract_metadata
  "$OSFM_ENV/bin/python" - "$D" <<'PY'
import json, sys
d = sys.argv[1]
cams = json.load(open(f"{d}/camera_models.json"))
out = {}
for k, c in cams.items():
    f = c.get("focal_x", c.get("focal"))
    out[k] = {"projection_type": "brown", "width": c["width"], "height": c["height"], "focal_x": f, "focal_y": f,
              "c_x": 0.0, "c_y": 0.0, "k1": 0.0, "k2": 0.0, "k3": 0.0, "p1": 0.0, "p2": 0.0}
json.dump(out, open(f"{d}/camera_models_overrides.json", "w"), indent=1)
print("cameras:", {k: round(v["focal_x"], 4) for k, v in out.items()})
PY
fi
[ -f "$D/reconstruction.json" ] || for s in extract_metadata detect_features match_features create_tracks reconstruct; do osfm $s; done
source "$HERE/env.sh"
$PY "$HERE/kjr/osfm_import.py" "$D" "$W/odm"
