#!/bin/bash
# 패널 영역 풀: 학습 정상 15 · 암 30.
#   원본은 pool15.sh · pool15_rand.sh · pool15_randb.sh 셋이었다. 세 파일의 차이는
#   「어느 패널 목록을 쓰는가」 하나뿐이어서 인자로 받는다. 값은 그대로 옮겼다.
#
#   실제 패널 풀과 귀무 풀은 자리를 가른다. 귀무를 만들다 실제 풀을 덮으면
#   학습을 처음부터 다시다. 칸마다 뽑히는 리드는 풀 구성과 무관하므로
#   가르더라도 비교 가능성은 유지된다.
#
#   사용
#     bash run_pool.sh cellline     # 실제 패널 (panel_union3_cellline.csv -> j_pool<GEN>_cellline)
#     bash run_pool.sh rand       # 귀무A    (panel_union_rand3.csv     -> j_pool<GEN>_rand)
#     bash run_pool.sh randb      # 귀무B    (panel_union_randb3.csv    -> j_pool<GEN>_randb)
set -eu

_HERE=$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)
# 환경변수로 주는 값이 이긴다. meth_config.py 의 값() 과 같은 순서로 맞춘다.
#   conf 를 그냥 source 하면 환경변수를 덮어써서, 환경변수로 격리한 실행이
#   조용히 실 경로로 넘어간다. 미리 보관하고 되돌린다.
#   값은 eval 에 태우지 않는다. 이름만 고정 목록으로 돈다.
# 보관 대상을 손으로 들고 있으면 목록에 없는 키는 conf 가 환경변수를
#   덮는다. meth_config.py 의 「환경변수가 이긴다」와 반대가 된다. conf 에 적힌 키를
#   「전부」 보관한다. 목록을 유지할 필요가 없어진다.
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
: "${METH_ROOT:=/ssd_data/Methylation}"
# 셸이 기본값을 채우면 conf 를 비웠을 때 관문이 못 잔다.
#   정본은 config.conf 와 meth_config.py 둘이다. 셸은 값을 만들지 않는다.
: "${GEN:?GEN 이 없습니다 — 0_setup/config.conf 의 GEN 을 채우십시오}"
: "${DMR_OUT:=$METH_ROOT/results/dmr}"
# use15_* 로 「세대 15」 를 박고 있었다. GEN=16 이면 16세대 풀을
#   15세대 목록으로 만든다. 세대가 섞이는데 아무 말도 안 난다.
: "${USE_NORMAL:=$HOME/tmp/use${GEN}_normal.txt}"   # 쓸 정상 BAM 목록
: "${USE_GBM:=$HOME/tmp/use${GEN}_gbm.txt}"         # 쓸 암 BAM 목록

WHICH=${1:-}
case "$WHICH" in
  cellline) PANEL=$DMR_OUT/j_panel$GEN/panel_union3_cellline.csv ;;
  rand)   PANEL=$DMR_OUT/j_panel$GEN/panel_union_rand3.csv ;;
  randb)  PANEL=$DMR_OUT/j_panel$GEN/panel_union_randb3.csv ;;
  *) echo "사용: $(basename "$0") cellline|rand|randb"; exit 1 ;;
esac
OUT=$DMR_OUT/j_pool${GEN}_$WHICH

[ -f "$PANEL" ] || { echo "패널 목록이 없다: $PANEL"; exit 1; }
for f in "$USE_NORMAL" "$USE_GBM"; do
  [ -f "$f" ] || { echo "쓸 BAM 목록이 없다: $f"; exit 1; }
done

say(){ echo "[$(date +%H:%M:%S)] $*"; }
say "패널 $PANEL"
say "산출 $OUT"

say "정상 -> pool_normal_all"
python -u "$_HERE/11_pool.py" --panel "$PANEL" --use "$USE_NORMAL" \
        --out "$OUT/pool_normal_all.bam"

say "암 -> pool_gbm_all"
python -u "$_HERE/11_pool.py" --panel "$PANEL" --use "$USE_GBM" \
        --out "$OUT/pool_gbm_all.bam"

say "완료"
