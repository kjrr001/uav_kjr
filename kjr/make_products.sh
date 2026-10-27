#!/usr/bin/env bash
# ODM-KJR last step: deliverables from <dir with ortho.tif + dsm_fused.tif> into <out dir>.
#   <B>_kjr_ortho.tif        COG RGBA DEFLATE
#   <B>_kjr_dsm.tif          COG Float32, nodata -9999, on the ortho grid (bilinear from the fused DSM)
#   <B>_kjr_{ortho,dsm}_envi.tif  uncompressed striped BigTIFF Float32 NaN (ENVI)
#   <B>_kjr_ortho.jpg, _ortho_85pct.jpg, _dsm_hillshade.jpg (+ .wld/.aux.xml)
# No DTM: the site has no vegetation, DSM = DTM.
# Usage: make_products.sh <BLOCK e.g. KJ1> <dir with ortho.tif, dsm_fused.tif> <out dir> [<name prefix>, default <BLOCK>_kjr]
set -euo pipefail
B="$1"; SRC="$2"; OUT="$3"
PYG="${PY:-python3}"   # any python with GDAL bindings
P="${4:-${B}_kjr}"
mkdir -p "$OUT"
TMP=$(mktemp -d -p "$OUT" .tmp.XXXX)
trap 'rm -rf "$TMP"' EXIT

echo "== ortho COG"
gdal_translate -q -of COG -co COMPRESS=DEFLATE -co BIGTIFF=YES "$SRC/ortho.tif" "$OUT/${P}_ortho.tif"

echo "== DSM on the ortho grid"
read -r XMIN YMAX XRES W H <<<"$($PYG -c "
from osgeo import gdal; d=gdal.Open('$SRC/ortho.tif'); g=d.GetGeoTransform()
print(g[0], g[3], g[1], d.RasterXSize, d.RasterYSize)")"
XMAX=$($PYG -c "print($XMIN + $XRES * $W)"); YMIN=$($PYG -c "print($YMAX - $XRES * $H)")
gdalwarp -q -r bilinear -te "$XMIN" "$YMIN" "$XMAX" "$YMAX" -ts "$W" "$H" \
  -srcnodata nan -dstnodata -9999 -ot Float32 -co TILED=YES -co BIGTIFF=YES \
  "$SRC/dsm_fused.tif" "$TMP/dsm.tif"
gdal_translate -q -of COG -co COMPRESS=DEFLATE -co PREDICTOR=YES -co BIGTIFF=YES "$TMP/dsm.tif" "$OUT/${P}_dsm.tif"

echo "== ENVI-ready rasters"
$PYG - "$OUT/${P}_ortho.tif" "$OUT/${P}_dsm.tif" "$OUT/${P}_ortho_envi.tif" "$OUT/${P}_dsm_envi.tif" <<'EOF'
import sys
import numpy as np
from osgeo import gdal
gdal.UseExceptions()
o_src, d_src, o_out, d_out = sys.argv[1:5]
STEP = 512
def create(src, path, nb):
    d = gdal.GetDriverByName("GTiff").Create(path, src.RasterXSize, src.RasterYSize, nb, gdal.GDT_Float32,
                                             ["COMPRESS=NONE", "TILED=NO", "BIGTIFF=YES"])
    d.SetGeoTransform(src.GetGeoTransform()); d.SetProjection(src.GetProjection())
    for i in range(nb): d.GetRasterBand(i + 1).SetNoDataValue(float("nan"))
    return d
s = gdal.Open(o_src); W, H = s.RasterXSize, s.RasterYSize; d = create(s, o_out, 3)
for y in range(0, H, STEP):
    n = min(STEP, H - y); m = s.GetRasterBand(4).ReadAsArray(0, y, W, n) == 0
    for i in range(3):
        a = s.GetRasterBand(i + 1).ReadAsArray(0, y, W, n).astype("float32"); a[m] = np.nan
        d.GetRasterBand(i + 1).WriteArray(a, 0, y)
d = None
s = gdal.Open(d_src); d = create(s, d_out, 1); b = s.GetRasterBand(1)
for y in range(0, H, STEP):
    n = min(STEP, H - y); a = b.ReadAsArray(0, y, W, n).astype("float32"); a[a <= -9000] = np.nan
    d.GetRasterBand(1).WriteArray(a, 0, y)
d.GetRasterBand(1).ComputeStatistics(False); d = None
print("envi written")
EOF

echo "== JPGs"
gdal_translate -q -of JPEG -b 1 -b 2 -b 3 -co QUALITY=92 -co WORLDFILE=YES "$OUT/${P}_ortho.tif" "$OUT/${P}_ortho.jpg"
gdal_translate -q -of JPEG -b 1 -b 2 -b 3 -outsize 85% 85% -co QUALITY=90 -co WORLDFILE=YES "$OUT/${P}_ortho.tif" "$OUT/${P}_ortho_85pct.jpg"
gdaldem hillshade -q -az 315 -alt 45 -compute_edges -co TILED=YES -co BIGTIFF=YES "$OUT/${P}_dsm.tif" "$TMP/hs.tif"
gdal_translate -q -of JPEG -co QUALITY=90 -co WORLDFILE=YES "$TMP/hs.tif" "$OUT/${P}_dsm_hillshade.jpg"
ls -la "$OUT"
