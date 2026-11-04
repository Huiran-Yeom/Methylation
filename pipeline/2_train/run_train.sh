#!/bin/bash
# 2단계: 학습
#
#   패널 영역에서 암, 정상 리드를 비율마다 섞고 피처 5개 만듦
#   피처 조합마다 모델 4개 교차검증해서 가장 좋은 것 저장
#
#   피처   mean, entropy, jsd, pdr, llr
#   비율   config.conf 의 MUT_RATIOS (0 ~ 5%)
#   모델   Logistic Regression, Random Forest, XGBoost, SVM
#
#   VERSION 하나씩 돎. 다른 판은 VERSION 바꿔서 다시 실행
#     VERSION=j15jsd200_5k_uni PANEL=jsd200 bash run_train.sh
#
#   사용법
#     bash run_train.sh            # 전 단계
#     bash run_train.sh 05         # 한 단계만
set -u

_HERE=$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)
S=$_HERE/steps
ONE=${1:-}

# config.conf 읽기. 환경변수가 우선 (meth_config.py 와 같은 규칙)
# conf 키를 미리 보관해 뒀다가 conf 읽은 뒤 되돌림
_CONFFILE=""
for _d in "$_HERE" "$_HERE/.." "$_HERE/../.." "$_HERE/../../.."; do
  [ -f "$_d/0_setup/config.conf" ] && { _CONFFILE="$_d/0_setup/config.conf"; break; }
done
_PREKEYS=$(awk -F= '/^[A-Za-z_][A-Za-z0-9_]*=/{print $1}' "${_CONFFILE:-/dev/null}" 2>/dev/null | tr '\n' ' ')
for _k in $_PREKEYS; do eval "_pre_$_k=\${$_k-}"; done
for _d in "$_HERE/.." "$_HERE/../.."; do
  if [ -f "$_d/0_setup/config.conf" ]; then
    set -a; . "$_d/0_setup/config.conf"; METH_CONF_DIR="$_d/0_setup"; set +a; break
  fi
done
for _k in $_PREKEYS; do
  eval "_v=\$_pre_$_k"
  if [ -n "${_v:-}" ]; then eval "export $_k=\$_pre_$_k"; fi
done
unset _v
: "${METH_CONF_DIR:?0_setup/config.conf 를 못 찾았다}"
# BLAS 스레드 고정 - 로지스틱 회귀 값이 스레드 수 따라 조금 달라짐
# 환경변수에 이미 있으면 그 값 사용, BLAS_THREADS=none 이면 고정 안 함
# 진입점(run_*.sh, make_report.sh)에서만 고정. steps/ 직접 실행하면 안 걸림
_bt=${BLAS_THREADS:-1}
case "$_bt" in
  [Nn][Oo][Nn][Ee]) ;;
  ''|*[!0-9]*|0) echo "멈춤: BLAS_THREADS 는 none 또는 1 이상의 정수여야 합니다: $_bt" >&2; exit 1 ;;
  *) export OMP_NUM_THREADS=${OMP_NUM_THREADS:-$_bt} OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-$_bt} MKL_NUM_THREADS=${MKL_NUM_THREADS:-$_bt}
     [ "$OMP_NUM_THREADS" = "$OPENBLAS_NUM_THREADS" ] && [ "$OMP_NUM_THREADS" = "$MKL_NUM_THREADS" ]        || echo "주의: 스레드 값이 갈린다 · OMP=$OMP_NUM_THREADS OPENBLAS=$OPENBLAS_NUM_THREADS MKL=$MKL_NUM_THREADS" >&2 ;;
esac

say(){ echo "[$(date +%H:%M:%S)] $*"; }
step(){  # step <번호> <제목> <파일> [인자...]
  # 파일 뒤 인자(09 의 --allcombos --pooled 등)까지 전부 넘김
  local no=$1 title=$2; shift 2
  if [ -n "$ONE" ] && [ "$ONE" != "$no" ]; then return 0; fi
  say "── $no  $title"
  ( cd "$S" && python -u "$@" ) || { say "!! 실패 $no"; exit 1; }
}

# GEN 은 기본값 없음. conf 에 없으면 중단
: "${GEN:?GEN 이 없습니다 — 0_setup/config.conf 의 GEN 을 채우십시오}"
# VERSION, PANEL 은 09 에서 처음 쓰지만 몇 시간 뒤라 여기서 먼저 확인
: "${VERSION:?VERSION 이 없습니다 — config.conf 에 넣거나 VERSION=j15bl200_5k_uni 로 주십시오}"
: "${PANEL:?PANEL 이 없습니다 — config.conf 에 넣거나 PANEL=bl200 으로 주십시오}"
: "${METH_ROOT:=/ssd_data/Methylation}"

# 06, 07 이 읽는 후보, 참조 경로 (안 넘기면 세대 없는 기본 경로 보고 중단)
# 3_validate/steps/run_one_sample.sh 와 같은 값이어야 함
export LLR_REF=${LLR_REF:-$METH_ROOT/results/dmr/j_llrref${GEN}/llrref_K${PATTERN_K:-3}_minr${MIN_READS_PER_WINDOW:-6}.parquet}
export LLR_CAND=${LLR_CAND:-$METH_ROOT/results/dmr/j_candidates${GEN}}
export JSD_CAND=${JSD_CAND:-$METH_ROOT/results/dmr/j_candidates${GEN}}
export JSD_JCOUNT=${JSD_JCOUNT:-$METH_ROOT/results/dmr/j_jsdcount${GEN}}

say "학습 시작 — VERSION=${VERSION:-?} · PANEL=${PANEL:-?} · 깊이 ${TARGET_DEPTH_READS:-?}"
say "  후보 $JSD_CAND · 참조 $LLR_REF"

step 01 "DMR 영역 리드 추출·정렬"   01_Preprocessing_for_ML_RawData_Mixing.py
step 02 "비율마다 섞어 BAM 쓰기"     02_Preprocessing_for_ML_RawData_Mixing.py
step 03 "Bismark 메틸화 추출"        03_Preprocessing_for_ML_RawData_Mixing.py
step 04 "커버리지 확인"              04_coverage_check.py
step 05 "mean · entropy 행렬"        05_Feature_Matrix.py
step 06 "jsd · pdr 행렬"             06_jsdfeat.py
step 07 "llr 행렬 (theta 별 · LLR_max)" 07_llrfeat.py
step 08 "조합별 모델 넷 교차검증"    08_ML_classifier.py
# 09 는 --allcombos --pooled 필수 - 빼면 조합 2개만 저장, 다른 세대 판까지 덮어씀
# (4단계가 쓰는 mean|entropy|jsd|llr|pdr 조합이 안 생김)
step 09 "조합마다 BEST 굳히기"       09_save_best_j.py --allcombos --pooled "$VERSION"

say "학습 끝: <뿌리>/results/panels/<앞머리>/<판>/step04_ML_classifier/"
