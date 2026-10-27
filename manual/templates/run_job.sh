#!/usr/bin/env bash
# ======================================================================================
#  uav_kjr 표준 실행 템플릿                                         manual/templates/run_job.sh
# ======================================================================================
#  프로젝트 폴더에 복사해서 쓴다 (new_block.sh 가 복사해 준다). 위의 설정만 고친다.
#    ./run_job.sh              배경 실행 (nohup). 로그: logs/run_<시각>.log
#    ./run_job.sh fg           앞에서 실행 (터미널을 닫으면 멈춘다)
#    ./run_job.sh check        검사만 (--dry-run): 실행할 단계·명령을 보여 주고 아무것도 안 함
#    ./run_job.sh status       돌고 있는지, 블록별 마지막 진행 상황
#    ./run_job.sh stop         돌고 있는 작업 멈추기 (다시 ./run_job.sh 하면 이어서 한다: EXISTING skip)
# ======================================================================================

# ---------------------------------------------------------------- 설정
KJR_REPO=/work/eridu_drone/kjr/odm-kjr     # uav_kjr 저장소 (kjr_run 이 있는 곳)
CFGS=( *.cfg )                             # 실행할 블록 cfg (차례로). 예: CFGS=( JJ1.cfg JJ3.cfg )
EXCLUDE=( common.cfg )                     # 블록이 아닌 cfg
START=                                     # 비우면 각 cfg 의 START.  예: START=BBA
END=                                       # 비우면 각 cfg 의 END.    예: END=DIAG
EXISTING=                                  # 비우면 각 cfg 의 EXISTING. skip | archive | delete

# ---------------------------------------------------------------- 이하 고치지 않는다
set -euo pipefail
cd "$(dirname "$0")"
B=(); for c in "${CFGS[@]}"; do [[ " ${EXCLUDE[*]} " == *" $c "* ]] || B+=("$c"); done
[ ${#B[@]} -gt 0 ] || { echo "블록 cfg 가 없다 (CFGS)"; exit 1; }
OPT=(); [ -n "$START" ] && OPT+=(--start "$START"); [ -n "$END" ] && OPT+=(--end "$END")
[ -n "$EXISTING" ] && OPT+=(--existing "$EXISTING")
RUN=("$KJR_REPO/kjr_run" "${B[@]}" "${OPT[@]}")
mkdir -p logs
running() { [ -f logs/run.pid ] && kill -0 "$(cat logs/run.pid)" 2>/dev/null; }
case "${1:-bg}" in
  check)  "${RUN[@]}" --dry-run ;;
  fg)     running && { echo "이미 실행 중 (pid $(cat logs/run.pid))"; exit 1; }
          "${RUN[@]}" 2>&1 | tee "logs/run_$(date +%Y%m%d_%H%M%S).log" ;;
  bg)     running && { echo "이미 실행 중 (pid $(cat logs/run.pid))"; exit 1; }
          "${RUN[@]}" --dry-run > /dev/null || { echo "검사 실패: ./run_job.sh check 로 확인"; exit 1; }
          L="logs/run_$(date +%Y%m%d_%H%M%S).log"
          setsid nohup "${RUN[@]}" > "$L" 2>&1 < /dev/null &   # 새 프로세스 그룹: stop 이 하위 작업까지 멈춘다
          echo $! > logs/run.pid; ln -sf "$(basename "$L")" logs/latest.log
          echo "시작: pid $!, 블록 ${B[*]}"; echo "진행 보기: tail -f $PWD/logs/latest.log   또는 ./run_job.sh status" ;;
  status) if running; then echo "실행 중 (pid $(cat logs/run.pid))"; else echo "실행 중 아님"; fi
          [ -e logs/latest.log ] && { grep "^== " logs/latest.log | tail -5; grep -A30 "FAILED" logs/latest.log | head -40 || true; } ;;
  stop)   running && { P=$(cat logs/run.pid); kill -TERM -- "-$P" 2>/dev/null || kill "$P"; echo "멈춤 (pid $P 와 하위 작업)"; } \
                  || echo "실행 중 아님" ;;
  *)      sed -n 5,9p "$0"; exit 1 ;;
esac
