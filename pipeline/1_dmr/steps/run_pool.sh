#!/bin/bash
# 패널 영역 풀 BAM 만들기 (정상, 암 각각)
#   실제 패널, 귀무A, 귀무B 풀은 패널 목록만 다르고 산출 폴더도 따로 씀
#   풀이 다르면 2단계 믹싱 리드가 달라짐 - 다른 풀의 특징을 열로 이어 붙이지 않음
#
#   사용
#     bash run_pool.sh cellline     # 실제 패널 (panel_union3_cellline.csv -> j_pool<GEN>_cellline)
#     bash run_pool.sh rand       # 귀무A    (panel_union_rand3.csv     -> j_pool<GEN>_rand)
#     bash run_pool.sh randb      # 귀무B    (panel_union_randb3.csv    -> j_pool<GEN>_randb)
set -eu

_HERE=$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)
# config.conf 읽기. 환경변수가 우선 (meth_config.py 와 같은 규칙)
# conf 키를 미리 보관해 뒀다가 conf 읽은 뒤 되돌림
_CONFFILE=""
for _d in "$_HERE" "$_HERE/.." "$_HERE/../.." "$_HERE/../../.."; do
  [ -f "$_d/0_setup/config.conf" ] && { _CONFFILE="$_d/0_setup/config.conf"; break; }
done
_PREKEYS=$(awk -F= '/^[A-Za-z_][A-Za-z0-9_]*=/{print $1}' "${_CONFFILE:-/dev/null}" 2>/dev/null | tr '\n' ' ')
for _k in $_PREKEYS; do eval "_pre_$_k=\${$_k-}"; done
for _d in "$_HERE" "$_HERE/.." "$_HERE/../.." "$_HERE/../../.."; do
  if [ -f "$_d/0_setup/config.conf" ]; then
    set -a; . "$_d/0_setup/config.conf"; METH_CONF_DIR="$_d/0_setup"; set +a; break
  fi
done
for _k in $_PREKEYS; do
  eval "_v=\$_pre_$_k"
  if [ -n "${_v:-}" ]; then eval "export $_k=\$_pre_$_k"; fi
done
unset _v
: "${METH_CONF_DIR:?config.conf 를 못 찾았다 — 트리 밖에서 돌리려면 METH_CONF_DIR 를 직접 주라}"
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
: "${METH_ROOT:=/ssd_data/Methylation}"
# GEN 은 기본값 없음. conf 에 없으면 중단
: "${GEN:?GEN 이 없습니다 — 0_setup/config.conf 의 GEN 을 채우십시오}"
: "${DMR_OUT:=$METH_ROOT/results/dmr}"
# 검체 목록 이름에 GEN 붙임 (세대 섞임 방지)
: "${USE_NORMAL:=$HOME/tmp/use${GEN}_normal.txt}"   # 쓸 정상 BAM 목록
: "${USE_GBM:=$HOME/tmp/use${GEN}_gbm.txt}"         # 쓸 암 BAM 목록

# 패널 목록 경로는 PANEL_CSV 에. PANEL 은 conf 의 판 이름이라 덮어쓰면 meth_config 가드가 멈춤
WHICH=${1:-}
case "$WHICH" in
  cellline) PANEL_CSV=$DMR_OUT/j_panel$GEN/panel_union3_cellline.csv
            # 새 이름 없으면 옛 트리 이름(panel_union3_세포주.csv) 읽음
            [ -f "$PANEL_CSV" ] || PANEL_CSV=$DMR_OUT/j_panel$GEN/panel_union3_세포주.csv ;;
  rand)   PANEL_CSV=$DMR_OUT/j_panel$GEN/panel_union_rand3.csv ;;
  randb)  PANEL_CSV=$DMR_OUT/j_panel$GEN/panel_union_randb3.csv ;;
  *) echo "사용: $(basename "$0") cellline|rand|randb"; exit 1 ;;
esac
OUT=$DMR_OUT/j_pool${GEN}_$WHICH

[ -f "$PANEL_CSV" ] || { echo "패널 목록이 없다: $PANEL_CSV"; exit 1; }
for f in "$USE_NORMAL" "$USE_GBM"; do
  [ -f "$f" ] || { echo "쓸 BAM 목록이 없다: $f"; exit 1; }
done

say(){ echo "[$(date +%H:%M:%S)] $*"; }
say "패널 $PANEL_CSV"
say "산출 $OUT"

say "정상 -> pool_normal_all"
python -u "$_HERE/11_pool.py" --panel "$PANEL_CSV" --use "$USE_NORMAL" \
        --out "$OUT/pool_normal_all.bam"

say "암 -> pool_gbm_all"
python -u "$_HERE/11_pool.py" --panel "$PANEL_CSV" --use "$USE_GBM" \
        --out "$OUT/pool_gbm_all.bam"

say "완료"
