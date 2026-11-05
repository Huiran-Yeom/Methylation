#!/bin/bash
# 패널 크기 N 고르기: 학습 자료 안에서 단위째 5접기 CV
#
#   접기마다 학습 접기로 1단계 → 크기마다 jsd, bl 패널 → 학습 → 빼 둔 단위 채점
#   접기 AUC 로 방식마다 1-SE 규칙으로 N 을 고름. 검증 코호트는 안 봄
#
#   단계
#     folds   01 접기 나누기 (원자료는 링크만)
#     dmr     02 접기마다 1단계, 크기별 패널, 풀, 코호트 설정
#     cells   03 칸(접기 × jsd/bl × N)마다 학습, 채점
#     pick    04 고르기
#
#   사용법 (pipeline/ 에서)
#     NORMAL_SET=… GBM_COV_DIR=… NORMAL_BAM_DIR=… GBM_BAM_DIR=… bash size_cv/run_size_cv.sh
#     bash size_cv/run_size_cv.sh dmr 3        # 한 단계만, 접기 3 만
#     bash size_cv/run_size_cv.sh cell 1 jsd200   # 칸 하나만
#     CV_DRY=1 bash size_cv/run_size_cv.sh pick   # 표 형식, 칸 검사만
#
#   설정 (환경변수, 기본값)
#     CV_ROOT=$METH_ROOT/size_cv   산출 뿌리 (원자료 밖)
#     CV_SIZES="100 200 300 500 1000 2000"   CV_NFOLD=5
#     CV_RATIO=0.005  CV_DEPTH=5000  CV_NREP=10   접기 AUC 의 비율, 깊이, 단위마다 반복
#     CV_LLR_HELDOUT=0   접기 안 학습 풀을 다시 거르지 않음
#     CV_TOPCAP=N        07 --topcap. N 이면 크기와 같게
#     CV_STRATA=         층 csv (sample,stratum). 없으면 정상 / 암
#     CV_NCELL=2  CV_NV=4  CV_MEMW=28  CV_DISK_MIN=150   동시 칸, 칸 안 검체 동시, 메모리 문(G), 디스크 하한(G)
#     동시 수는 곱으로 셈: CV_NCELL × (NSLICE 또는 CV_NV). 기본 2 × 4 = 8
set -u

_HERE=$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)
P=$(cd "$_HERE/.." && pwd)
STEPS=$_HERE/steps
ONE=${1:-}; shift || true

# config.conf 읽기. 환경변수가 우선 (run_dmr.sh 와 같은 규칙)
_CONFFILE=$P/0_setup/config.conf
[ -f "$_CONFFILE" ] || { echo "멈춤: $_CONFFILE 없음" >&2; exit 1; }
_PREKEYS=$(awk -F= '/^[A-Za-z_][A-Za-z0-9_]*=/{print $1}' "$_CONFFILE" | tr '\n' ' ')
for _k in $_PREKEYS; do eval "_pre_$_k=\${$_k-}"; done
set -a; . "$_CONFFILE"; METH_CONF_DIR=$P/0_setup; set +a
for _k in $_PREKEYS; do
  eval "_v=\$_pre_$_k"
  if [ -n "${_v:-}" ]; then eval "export $_k=\$_pre_$_k"; fi
done
unset _v
: "${GEN:?GEN 이 없습니다}" "${METH_ROOT:?}"

export P GEN
export CV_ROOT=${CV_ROOT:-$METH_ROOT/size_cv}
export CV_SIZES=${CV_SIZES:-100 200 300 500 1000 2000} CV_NFOLD=${CV_NFOLD:-5}
export CV_RATIO=${CV_RATIO:-0.005} CV_DEPTH=${CV_DEPTH:-5000} CV_NREP=${CV_NREP:-10}
export CV_LLR_HELDOUT=${CV_LLR_HELDOUT:-0} CV_TOPCAP=${CV_TOPCAP:-N} CV_STRATA=${CV_STRATA:-}
export CV_NCELL=${CV_NCELL:-2} CV_NV=${CV_NV:-4} CV_MEMW=${CV_MEMW:-28} CV_DISK_MIN=${CV_DISK_MIN:-150}
case "$CV_TOPCAP" in N|[1-9]*) ;; *) echo "멈춤: CV_TOPCAP 은 N 또는 숫자: $CV_TOPCAP" >&2; exit 1 ;; esac
case "$(cd "$(dirname "$CV_ROOT")" 2>/dev/null && pwd -P)/" in
  "$(cd "$METH_ROOT" 2>/dev/null && pwd -P)/"*) ;;
  *) echo "멈춤: CV_ROOT 는 METH_ROOT 밑이어야 함 (정리 범위 제한): $CV_ROOT" >&2; exit 1 ;;
esac
# 설정 서명: conf, meth_config.py, 지금 값 (이 모듈이 정하는 키 빼고), 자료 위치. 02, 03 완료 표지에 들어감
_own=" METH_ROOT WORK_ROOT TMP_ROOT LLR_HELDOUT PANELS PANEL VERSION TARGET_DEPTH_READS COHORT COHORT_CONF VAL_COHORTS "
CV_CONF_SIG=$( { md5sum < "$_CONFFILE"; md5sum < "$P/0_setup/meth_config.py"
  for _k in $_PREKEYS; do case "$_own" in *" $_k "*) ;; *) echo "$_k=${!_k-}" ;; esac; done
  for _k in SID_RE GROUP_RE CANCER_LABEL NORMAL_SET GBM_COV_DIR NORMAL_BAM_DIR GBM_BAM_DIR; do echo "$_k=${!_k-}"; done
  } | md5sum | cut -c1-12)
export CV_CONF_SIG
mkdir -p "$CV_ROOT"
say(){ echo "[$(date +%F' '%T)] $*"; }

# 한 번만 돌게
exec 9> "$CV_ROOT/.lock"
flock -n 9 || { say "멈춤: 다른 실행이 $CV_ROOT 를 쓰는 중"; exit 1; }

do_folds(){
  python -u "$STEPS/01_folds.py" --out "$CV_ROOT" --nfold "$CV_NFOLD" ${CV_STRATA:+--strata "$CV_STRATA"} ${CV_DRY:+--dry}
}
do_dmr(){
  local ks=${*:-$(seq 1 "$CV_NFOLD")} k
  for k in $ks; do
    bash "$STEPS/02_fold_dmr.sh" "$k" > "$CV_ROOT/fold$k.dmr.log" 2>&1 || { say "멈춤: 접기 $k 1단계 실패 · $CV_ROOT/fold$k.dmr.log"; return 1; }
    say "접기 $k 1단계 끝"
  done
}
do_cells(){
  local ks=${*:-$(seq 1 "$CV_NFOLD")} k N m fail
  for k in $ks; do for N in $CV_SIZES; do for m in jsd bl; do echo "$k $m$N"; done; done; done \
    | xargs -P "$CV_NCELL" -L 1 bash -c 'bash "'"$STEPS"'/03_cells.sh" "$0" "$1" >> "'"$CV_ROOT"'/cells.log" 2>&1 || echo "f$0 $1"' \
    > "$CV_ROOT/cells_failed.txt"
  fail=$(grep -c . "$CV_ROOT/cells_failed.txt")
  [ "$fail" = 0 ] || { say "멈춤: 실패한 칸 $fail · $CV_ROOT/cells_failed.txt · 로그 $CV_ROOT/cells.log"; return 1; }
  say "칸 끝 · 계산 불가 $(ls "$CV_ROOT"/cells/*.infeasible 2>/dev/null | wc -l)"
}
do_pick(){
  python -u "$STEPS/04_pick.py" --root "$CV_ROOT" --gen "$GEN" --sizes "$CV_SIZES" --nfold "$CV_NFOLD" \
    --ratio "$CV_RATIO" --nrep "$CV_NREP" --depth "$CV_DEPTH" ${CV_DRY:+--dry} | tee "$CV_ROOT/pick${CV_DRY:+_dry}.txt"
  return "${PIPESTATUS[0]}"
}

say "size_cv · 뿌리 $CV_ROOT · 크기 $CV_SIZES · 접기 $CV_NFOLD · 깊이 $CV_DEPTH · 비율 $CV_RATIO · TOPCAP $CV_TOPCAP"
case "$ONE" in
  folds) do_folds ;;
  dmr)   do_dmr "$@" ;;
  cells) do_cells "$@" ;;
  cell)  bash "$STEPS/03_cells.sh" "$@" ;;
  pick)  do_pick ;;
  '')    do_folds && do_dmr && do_cells && do_pick ;;
  *)     echo "사용: bash run_size_cv.sh [folds|dmr|cells|cell|pick] [접기…]" >&2; exit 1 ;;
esac
