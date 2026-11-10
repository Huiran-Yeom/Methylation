#!/bin/bash
# size_cv 02: 접기 하나의 1단계
#   학습 접기 cov 로 1단계 01, 02, 03, 04, 05 (04b, 04c 는 jsdb 용이라 안 씀)
#   크기마다 06 (jsd), 07 (bl) → 02b 로 검증 형식 + 합집합
#   학습 풀 하나, 빼 둔 정상 풀 하나, 빼 둔 암 단위마다 풀 하나 → 빼 둔 단위 코호트 설정
# 사용법: bash 02_fold_dmr.sh <접기 번호>   (run_size_cv.sh 가 환경을 정해서 부름)
set -u
K=${1:?접기 번호}
: "${CV_ROOT:?}" "${CV_SIZES:?}" "${CV_TOPCAP:?}" "${CV_RATIO:?}" "${CV_NREP:?}" "${CV_CONF_SIG:?}" "${GEN:?}" "${P:?}"
F=$CV_ROOT/fold$K
[ -d "$F/train_normal_cov" ] && [ -d "$F/train_cancer_cov" ] || { echo "멈춤: 01 을 먼저 ($F)"; exit 1; }
say(){ echo "[$(date +%F' '%T)] f$K $*"; }
run(){ say "\$ $*"; "$@" || { say "실패 rc=$? — 중단"; exit 1; }; }

# 격리: 바깥 값이 끼어들지 않게
unset JSD_CAND LLR_CAND USE_NORMAL USE_GBM USE_CANCER POOL_DIR COLSET_FROM SCORE_UNIT CONVERGENCE_GATE
export METH_ROOT=$F/root WORK_ROOT=$F/work TMP_ROOT=$F/tmp LLR_HELDOUT=$CV_LLR_HELDOUT
export NORMAL_SET=$F/train_normal_cov CANCER_COV_DIR=$F/train_cancer_cov GBM_COV_DIR=$F/train_cancer_cov   # 새 이름, 옛 이름 둘 다
export PANELS=$(for N in $CV_SIZES; do printf 'jsd%s,bl%s,' "$N" "$N"; done | sed 's/,$//')
mkdir -p "$WORK_ROOT" "$TMP_ROOT"
D=$METH_ROOT/results/dmr

# 이미 같은 코드, 같은 설정으로 끝난 접기면 건너뜀
SIG=$( { cat "$P"/size_cv/steps/02_fold_dmr.sh "$P"/size_cv/steps/02b_export_sizes.py "$P"/1_dmr/run_dmr.sh \
          "$P"/1_dmr/steps/0[1-7]*.py "$P"/1_dmr/steps/11_pool.py \
          "$P"/1_dmr/steps/moderated_t.py "$P"/1_dmr/steps/sample_id.py "$P"/1_dmr/steps/dmr_common.py "$P"/1_dmr/steps/jsd_resample.py \
          "$CV_ROOT/folds.csv" "$F"/heldout_*.txt | md5sum | cut -c1-12
        echo "$CV_SIZES|$CV_TOPCAP|$CV_RATIO|$CV_NREP|$CV_LLR_HELDOUT|$GEN|$CV_CONF_SIG"; } | tr '\n' ' ')
if [ -s "$F/.dmr_done" ] && [ "$(cat "$F/.dmr_done")" = "$SIG" ]; then say "이미 끝남 — 건너뜀"; exit 0; fi
if [ -e "$METH_ROOT" ] && [ -n "$(ls -A "$METH_ROOT" 2>/dev/null)" ]; then
  say "멈춤: $METH_ROOT 가 비어있지 않은데 완료 표지가 없거나 다름. 이 모듈이 만든 폴더이니 지우고 다시 돌리거나 그대로 쓸지 확인"
  exit 1
fi
say "시작 · 정상 $(ls "$NORMAL_SET"/*.cov.gz | wc -l) · 암 $(ls "$CANCER_COV_DIR"/*.cov.gz | wc -l) · 크기 $CV_SIZES · TOPCAP $CV_TOPCAP"

for s in 01 02 03 04 05; do run bash "$P/1_dmr/run_dmr.sh" $s; done
JF=$(ls "$D/j_jsd$GEN"/jsd_K${PATTERN_K:-3}_minr${MIN_READS_PER_WINDOW:-6}.parquet)
for N in $CV_SIZES; do
  T=$N; [ "$CV_TOPCAP" = N ] || T=$CV_TOPCAP
  run python -u "$P/1_dmr/steps/06_select.py" --cand "$D/j_candidates$GEN" --jsd "$JF" --out "$D/j_panel$GEN" --n "$N"
  run python -u "$P/1_dmr/steps/07_panel.py" --cand "$D/j_candidates$GEN" --out "$D/j_panel_dmr$GEN/bl${N}_cap$T" \
      --topn "$N" --topcap "$T"
done
run python -u "$P/size_cv/steps/02b_export_sizes.py" --dmr "$D" --gen "$GEN" --sizes "$CV_SIZES" --topcap "$CV_TOPCAP"

# 학습 풀: 학습 접기 정상, 암 전부 (LLR_HELDOUT=0 이면 11_pool 이 다시 거르지 않음)
S=$D/j_candidates$GEN/samples.csv
awk -F, 'NR>1 && $3=="Normal"{print $1}' "$S" > "$F/use_normal.txt"
awk -F, 'NR>1 && $3!="Normal"{print $1}' "$S" > "$F/use_cancer.txt"
say "학습 풀 목록: 정상 $(wc -l < "$F/use_normal.txt") · 암 $(wc -l < "$F/use_cancer.txt")"
PL=$D/j_pool${GEN}_cellline
run python -u "$P/1_dmr/steps/11_pool.py" --panel "$D/union.csv" --use "$F/use_normal.txt" --out "$PL/pool_normal_all.bam"
run python -u "$P/1_dmr/steps/11_pool.py" --panel "$D/union.csv" --use "$F/use_cancer.txt" --out "$PL/pool_gbm_all.bam"

# 빼 둔 단위: 정상 풀 하나, 암 단위마다 풀 하나, 코호트 설정
#   빼 둔 검체는 학습 samples.csv 에 없으니 11_pool 의 군 거르기를 끔 (LLR_HELDOUT=0)
VI=$METH_ROOT/val_in; mkdir -p "$VI"
run env LLR_HELDOUT=0 python -u "$P/1_dmr/steps/11_pool.py" --panel "$D/union.csv" --use "$F/heldout_normal.txt" --out "$VI/heldout_normal_pool.bam"
while read -r u; do
  [ -n "$u" ] || continue
  awk -F, -v u="$u" 'NR>1 && $2==u && $3!="Normal"{print $1}' "$CV_ROOT/folds.csv" > "$VI/use_$u.txt"
  run env LLR_HELDOUT=0 python -u "$P/1_dmr/steps/11_pool.py" --panel "$D/union.csv" --use "$VI/use_$u.txt" --out "$VI/${u}_pool.bam"
done < "$F/heldout_cancer_units.txt"
{
  echo "# size_cv 접기 $K 빼 둔 단위 코호트 (이름이 든 파일, CV_ROOT 밖으로 내지 않음)"
  echo "NORMAL=$VI/heldout_normal_pool.bam"
  echo "CANCER=$VI/{sample}_pool.bam"
  echo "SAMPLES=$(grep . "$F/heldout_cancer_units.txt" | tr '\n' ' ' | sed 's/ $//')"
  echo "NREP=$CV_NREP"
  echo "RATIOS=0.0, $CV_RATIO"
  echo "LAYOUT=mixed"
} > "$F/cohort.conf"
say "코호트 설정: 암 단위 $(grep -c . "$F/heldout_cancer_units.txt")"
printf '%s' "$SIG" > "$F/.dmr_done.tmp" && mv "$F/.dmr_done.tmp" "$F/.dmr_done"
say "끝"
