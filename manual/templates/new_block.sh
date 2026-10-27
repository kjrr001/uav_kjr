#!/usr/bin/env bash
# ======================================================================================
#  새 블록 시작: 프로젝트 폴더에 표준 템플릿으로 cfg 를 만든다       manual/templates/new_block.sh
# ======================================================================================
#  사용: new_block.sh <프로젝트 폴더> <BLOCK> <사진 폴더> <GCP 파일> ["설명"]
#   - 프로젝트 폴더가 없으면 만든다. common.cfg, run_job.sh 가 없으면 템플릿을 복사한다.
#   - <BLOCK>.cfg 를 block.cfg 템플릿으로 만든다 (이미 있으면 덮어쓰지 않는다).
#   - 사진 수, GCP 파일 형식, cfg 를 검사한다 (kjr_run --dry-run).
#  작업·제품 폴더 기본값: <프로젝트 폴더>/work/<BLOCK>, <프로젝트 폴더>/products/<BLOCK>
# ======================================================================================
set -euo pipefail
[ $# -ge 4 ] || { sed -n 5,9p "$0"; exit 1; }
T=$(cd "$(dirname "$0")" && pwd); REPO=$(cd "$T/../.." && pwd)
P=$(realpath -m "$1"); B="$2"; IMG=$(realpath "$3"); GCP=$(realpath "$4"); TITLE="${5:-$B}"
[[ "$B" =~ ^[A-Za-z0-9_-]+$ ]] || { echo "BLOCK 은 영문·숫자·_- 만: $B"; exit 1; }
N=$(find "$IMG" -maxdepth 1 -type f \( -iname '*.jpg' -o -iname '*.jpeg' \) | wc -l)
[ "$N" -gt 0 ] || { echo "사진(*.jpg)이 없다: $IMG"; exit 1; }
H=$(head -1 "$GCP"); G=$(awk 'NR>1 && NF>=6' "$GCP" | wc -l); GP=$(awk 'NR>1 && NF>=7 {print $7}' "$GCP" | sort -u | wc -l)
[[ "$H" =~ ^(EPSG:|WGS84|\+proj) ]] || echo "주의: GCP 첫 줄이 좌표계가 아닌 것 같다: $H"
mkdir -p "$P"
for f in common.cfg run_job.sh; do [ -e "$P/$f" ] || { cp "$T/$f" "$P/"; echo "복사: $P/$f"; }; done
sed -i "s|^KJR_REPO=.*|KJR_REPO=$REPO     # uav_kjr 저장소 (kjr_run 이 있는 곳)|" "$P/run_job.sh"
C="$P/$B.cfg"
if [ -e "$C" ]; then echo "이미 있음 (그대로 둠): $C"; else
  sed -e "s|<설명, 예: Eridu JJ-3 2026-07 촬영>|$TITLE|" -e "s|<BLOCK>|$B|" \
      -e "s|<INPUT_DIR>|$IMG|" -e "s|<GCP_FILE>|$GCP|" "$T/block.cfg" > "$C"
  echo "만듦: $C"
fi
echo "사진 $N 장, GCP 관측 $G 줄 (점 $GP 개), 좌표계 $H"
"$REPO/kjr_run" "$C" --dry-run > /tmp/new_block_$$.txt 2>&1 && echo "검사 통과: $(grep ' start' /tmp/new_block_$$.txt | sed 's/.*| stages/stages/' | cut -d'|' -f1)" \
  || { cat /tmp/new_block_$$.txt; rm -f /tmp/new_block_$$.txt; exit 1; }
rm -f /tmp/new_block_$$.txt
echo "다음: cd $P && ./run_job.sh check && ./run_job.sh"
