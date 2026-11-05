#!/bin/bash
# size_cv 03: 칸 하나 = 접기 × 방식(jsd, bl) × 크기 N
#   학습 run/lib.sh train_cell → 빼 둔 단위 채점 run_validate.sh
#   학습 완료 표지가 같은 코드, 같은 설정, 같은 02 결과면 학습만 건너뜀
#   표지가 없거나 다르면 이 칸의 옛 학습, 채점 산출을 지우고 새로 함 (2단계 02, 03 은 있는 파일을 이어 씀)
#   채점은 셸에서 안 건너뜀 (score_one_sample 이 모델 md5 로 판단)
# 사용법: bash 03_cells.sh <접기> <패널 이름 예: jsd200>   (run_size_cv.sh 가 환경을 정해서 부름)
set -u
K=${1:?접기}; PN=${2:?패널}
: "${CV_ROOT:?}" "${CV_DEPTH:?}" "${CV_NV:?}" "${CV_MEMW:?}" "${CV_DISK_MIN:?}" "${CV_CONF_SIG:?}" "${GEN:?}" "${P:?}"
F=$CV_ROOT/fold$K
say(){ echo "[$(date +%F' '%T)] f$K $PN $*"; }
[ -s "$F/.dmr_done" ] || { say "멈춤: 02 를 먼저"; exit 1; }

unset JSD_CAND LLR_CAND USE_NORMAL USE_GBM USE_CANCER POOL_DIR COLSET_FROM SCORE_UNIT CONVERGENCE_GATE
case "$CV_DEPTH" in *000) DG=$((CV_DEPTH / 1000))k ;; *) DG=$CV_DEPTH ;; esac
V=j${GEN}${PN}_${DG}_uni
export METH_ROOT=$F/root WORK_ROOT=$F/work TMP_ROOT=$F/tmp LLR_HELDOUT=$CV_LLR_HELDOUT
export PANELS=$PN
export COHORT=cv_f$K COHORT_CONF=$F/cohort.conf VAL_COHORTS=cv_f$K
export PANEL=$PN VERSION=$V TARGET_DEPTH_READS=$CV_DEPTH      # 채점용. train_cell 도 같은 값으로 다시 정함
mkdir -p "$WORK_ROOT" "$TMP_ROOT" "$CV_ROOT/cells"
PC=$METH_ROOT/results/dmr/j_panel_dmr${GEN}_panels/$PN/DMR_confirmed_$PN.csv
[ -s "$PC" ] && [ -s "$F/cohort.conf" ] || { say "멈춤: 패널 또는 코호트 설정 없음"; exit 1; }
MK=$CV_ROOT/cells/f${K}_$PN

# 학습: run/lib.sh 의 train_cell (01, 02 → 03 조각 NSLICE → 04~09 → BEST 31 확인)
_path=$PATH; . "$P/run/lib.sh"; PATH=$_path      # lib.sh 가 붙이는 서버 conda 경로는 안 씀
SIG=$( { cat "$P"/run/lib.sh "$P"/2_train/run_train.sh "$P"/2_train/steps/*.py "$PC" | md5sum | cut -c1-12
        echo "$CV_DEPTH|$CV_LLR_HELDOUT|$GEN|$CV_CONF_SIG|$(cat "$F/.dmr_done")"; } | tr '\n' ' ')
VD=$METH_ROOT/results/panels/j${GEN}${PN}/$V
if [ -s "$MK.train_done" ] && [ "$(cat "$MK.train_done")" = "$SIG" ] && [ "$(best_n "$VD")" = "$EXPECT_BEST" ]; then
  say "학습 이미 끝남 — 건너뜀"
elif [ -s "$MK.infeasible" ] && [ "$(cat "$MK.infeasible")" = "$SIG" ]; then
  say "계산 불가(풀 부족) 표지 있음 — 건너뜀"; exit 0
else
  RP=$(cd "$CV_ROOT" && pwd -P)
  for d in "$VD" "$METH_ROOT/sample_data/panels/j${GEN}${PN}/$V" "$METH_ROOT"/cv_f${K}_val/*/results/"$V" \
           "$METH_ROOT"/cv_f${K}_val/*/sample_data/panels/*/"$V" "$WORK_ROOT/cv_f${K}_val/$V"; do
    [ -e "$d" ] || continue
    case "$(cd "$d" && pwd -P)" in "$RP"/*) rm -rf -- "$d" ;; *) say "멈춤: CV_ROOT 밖 $d"; exit 1 ;; esac
  done
  rm -f "$MK.train_done" "$MK.infeasible"
  ( flock 8; while [ "$(free -g | awk 'NR==2{print $3}')" -ge "$CV_MEMW" ]; do sleep 60; done; sleep "${CV_STAGGER:-180}" ) \
    8> "$CV_ROOT/cells/start.lock"
  L=$TMP_ROOT/train_${PN}_$DG.log
  say "학습 시작 · 패널 $(($(wc -l < "$PC") - 1))블록 · $V"
  train_cell "$PN" "$DG" "$CV_DEPTH"; rc=$?
  if [ "$rc" != 0 ]; then
    if grep -q "풀 부족으로 중단" "$L"; then
      printf '%s' "$SIG" > "$MK.infeasible.tmp" && mv "$MK.infeasible.tmp" "$MK.infeasible"
      say "계산 불가(풀 부족)"; exit 0
    fi
    say "학습 실패 rc=$rc · $L"; exit 1
  fi
  printf '%s' "$SIG" > "$MK.train_done.tmp" && mv "$MK.train_done.tmp" "$MK.train_done"
  rm -rf -- "$METH_ROOT/sample_data/panels/j${GEN}${PN}/$V"     # 학습 섞기 중간물 (이 모듈이 만든 것)
  say "학습 끝"
fi

# 채점
VL=$TMP_ROOT/val_cv_f${K}_$V.log
say "채점 시작 · 검체 동시 $CV_NV"
bash "$P/3_validate/run_validate.sh" "$V" "$CV_NV" > "$VL" 2>&1; rc=$?
say "채점 rc=$rc · $(grep -o '채점된 검체 [0-9]*/[0-9]*' "$VL" | tail -1)"
[ "$rc" = 0 ] || exit 1

# 디스크가 모자라면 채점 끝난 칸의 섞은 BAM 정리 (CV_ROOT 밑, 이 모듈이 만든 것만)
av=$(df -BG --output=avail "$CV_ROOT" | tail -1 | tr -dc 0-9)
if [ "${av:-0}" -lt "$CV_DISK_MIN" ]; then
  for d in "$METH_ROOT"/cv_f${K}_val/*/results/"$V"/step03_Preprocessing_for_ML; do
    case "$(cd "$d" 2>/dev/null && pwd -P)" in "$(cd "$CV_ROOT" && pwd -P)"/*) rm -rf -- "$d" ;; esac
  done
  say "디스크 ${av}G — 섞은 BAM 정리"
fi
say "끝"
