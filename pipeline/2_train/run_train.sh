#!/bin/bash
# ── 2단계 · 학습 ──────────────────────────────────────────────────────
#
#   패널 영역에서 암·정상 조각을 비율마다 섞고(in-silico mixing), 피처 다섯을
#   만들고, 조합마다 모델 넷을 교차검증해 가장 좋은 것을 굳힌다.
#
#   피처 다섯   mean · entropy · jsd · pdr · llr
#   비율        config.conf 의 MUT_RATIOS (0 ~ 5%)
#   모델        Logistic Regression · Random Forest · XGBoost · SVM
#
#   판 하나(VERSION)에 대해 돈다. 여러 판을 돌리려면 VERSION 을 바꿔 다시 부른다:
#     VERSION=j15jsd200_5k_uni PANEL=jsd200 bash run_train.sh
#
#   사용
#     bash run_train.sh            # 전 단계
#     bash run_train.sh 05         # 그 단계만
set -u

_HERE=$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)
S=$_HERE/steps
ONE=${1:-}

# 2026-09-28: 보관 대상을 손으로 들고 있으면 목록에 없는 키는 conf 가 환경변수를
#   덮는다. meth_config.py 의 「환경변수가 이긴다」와 반대가 된다. conf 에 적힌 키를
#   「전부」 보관한다. 목록을 유지할 필요가 없어진다.
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

say(){ echo "[$(date +%H:%M:%S)] $*"; }
step(){  # step <번호> <제목> <파일> [인자...]
  # 2026-09-29: $3 만 쓰면 뒤에 붙인 플래그가 「조용히 버려진다」.
  #   09 에 --allcombos --pooled 를 붙였는데 그대로 사라졌다. 전부 넘긴다.
  local no=$1 title=$2; shift 2
  if [ -n "$ONE" ] && [ "$ONE" != "$no" ]; then return 0; fi
  say "── $no  $title"
  ( cd "$S" && python -u "$@" ) || { say "!! 실패 $no"; exit 1; }
}

# 2026-09-29: 아래 export 가 GEN·METH_ROOT 를 처음 읽는 자리다. 없으면 set -u 로
#   맨 「unbound variable」 이 난다. 다른 드라이버처럼 뜻이 있는 말로 멈춘다.
: "${GEN:?GEN 이 없습니다 — 0_setup/config.conf 의 GEN 을 채우십시오}"
: "${METH_ROOT:=/ssd_data/Methylation}"

# 2026-09-29: 06·07 이 보는 후보·참조 경로를 「여기서」 내보낸다.
#   예전에는 아무것도 안 내보내서 둘이 자기 기본값(세대 없는 j_candidates 등)으로
#   내려갔고, 그 둘의 세대 관문이 즉시 멈췄다. 문서대로 돌리면 06 에서 끝났다.
#   검증 러너(3_validate/steps/run_one_sample.sh)와 「같은 값」 이어야 한다.
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
# 2026-09-29: 인자 없이 부르면 09 는 조합 「둘」(mean|entropy, +jsd)만 굳히고,
#   판 목록을 results/panels/* 에서 「모든 세대에 걸쳐」 훑어 덮어쓴다.
#   그러면 사전등록 1차 조합 mean|entropy|jsd|llr|pdr 이 아예 안 만들어져
#   4단계가 「조합이 없음」 으로 죽는다. 실제 실행은 --allcombos --pooled 였다
#   (서버 산출물에 조합 31개·pooled 만 있는 것으로 확인).
step 09 "조합마다 BEST 굳히기"       09_save_best_j.py --allcombos --pooled "$VERSION"

say "학습 끝: <뿌리>/results/panels/<앞머리>/<판>/step04_ML_classifier/"
