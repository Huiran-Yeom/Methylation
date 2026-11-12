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
# 공통 머리: config.conf 읽기 (환경변수 우선), BLAS 스레드 고정, 암 쪽 자료 이름 맞추기 (0_setup/load_config.sh)
_LC=""
for _d in "$_HERE" "$_HERE/.." "$_HERE/../.." "$_HERE/../../.."; do
  [ -f "$_d/0_setup/load_config.sh" ] && { _LC="$_d/0_setup/load_config.sh"; break; }
done
: "${_LC:=${METH_CONF_DIR:+$METH_CONF_DIR/load_config.sh}}"
[ -f "${_LC:-}" ] || { echo "멈춤: 0_setup/load_config.sh 를 못 찾았다. 트리 밖에서 돌리려면 METH_CONF_DIR 를 직접 주라" >&2; exit 1; }
. "$_LC"
: "${METH_ROOT:=/ssd_data/Methylation}"
# GEN 은 기본값 없음. conf 에 없으면 중단
: "${GEN:?GEN 이 없습니다. 0_setup/config.conf 의 GEN 을 채우십시오}"
: "${DMR_OUT:=$METH_ROOT/results/dmr}"
# 검체 목록 이름에 GEN 붙임 (세대 섞임 방지)
: "${USE_NORMAL:=$DMR_OUT/j_panel$GEN/use${GEN}_normal.txt}"   # 쓸 정상 BAM 목록 (10c 가 씀, run_dmr 과 같은 자리)
: "${USE_GBM:=$DMR_OUT/j_panel$GEN/use${GEN}_gbm.txt}"         # 쓸 암 BAM 목록

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
