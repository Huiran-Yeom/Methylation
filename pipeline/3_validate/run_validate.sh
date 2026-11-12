#!/bin/bash
# 3단계: 외부검증
#
#   학습에 안 쓴 정상 배경에 외부 코호트 암 검체를 학습과 같은 비율, 깊이, 섞기로 섞고
#   학습된 모델로 채점. 바뀌는 건 입력 BAM 과 출력 경로뿐
#
#   학습 스크립트를 검체 폴더로 복사해서 돌림
#   2_train/steps 파일 이름이 config.conf 의 STEP_MAP 과 묶여 있음 (이름 바꾸면 같이 수정)
#
#   흐름
#     검체마다   make_val_config -> run_one_sample (피처 5종) -> score_one_sample (채점)
#     끝난 뒤    4_report/make_report.sh 가 검체별 채점 파일 읽어 표 만듦
#
#   사용법
#     bash run_validate.sh <판이름> [동시실행수]
#     COHORT=<your_cohort> bash run_validate.sh j15bl200_5k_uni 8
set -u

_HERE=$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)
S=$_HERE/steps

# 공통 머리: config.conf 읽기 (환경변수 우선), BLAS 스레드 고정, 암 쪽 자료 이름 맞추기 (0_setup/load_config.sh)
_LC=""
for _d in "$_HERE/.." "$_HERE/../.."; do
  [ -f "$_d/0_setup/load_config.sh" ] && { _LC="$_d/0_setup/load_config.sh"; break; }
done
: "${_LC:=${METH_CONF_DIR:+$METH_CONF_DIR/load_config.sh}}"
[ -f "${_LC:-}" ] || { echo "멈춤: 0_setup/load_config.sh 를 못 찾았다. 트리 밖에서 돌리려면 METH_CONF_DIR 를 직접 주라" >&2; exit 1; }
. "$_LC"
: "${COHORT:?COHORT 를 주십시오 (예: COHORT=<your_cohort>)}"

# flock 없으면 중단 - 없으면 전부 다른 실행이 잡은 걸로 잘못 셈 (Git Bash 에는 없음)
command -v flock >/dev/null 2>&1 || {
  echo "멈춤: flock 이 없습니다 (util-linux). 이 파이프라인은 검체별 잠금에 씁니다." >&2
  echo "      없으면 「빠짐」과 「남이 쥐고 있음」을 구분하지 못합니다." >&2
  exit 1; }

# COHORT 가 VAL_COHORTS 에 없으면 중단 (오타, 옛 코호트 방지)
# VAL_COHORTS 비어 있으면 아무 코호트나 허용. 쉼표 뒤 공백은 지움
if [ -n "${VAL_COHORTS:-}" ]; then
  _VC=$(printf %s "$VAL_COHORTS" | tr -d " 	")
  case ",${_VC}," in
    *",${COHORT},"*) : ;;
    *) echo "멈춤: COHORT=$COHORT 는 VAL_COHORTS 에 없습니다 ($VAL_COHORTS)" >&2
       echo "   오타이거나 옛 코호트입니다. 쓰려면 0_setup/config.conf 의" >&2
       echo "   VAL_COHORTS 에 더하십시오: 코드를 고치지 마십시오." >&2
       exit 1 ;;
  esac
fi

V=${1:?사용: bash run_validate.sh <판이름> [동시실행수]}
PAR=${2:-8}
# 학습 스크립트 원본 자리, make_val_config 자리
export TRAIN_SRC=${TRAIN_SRC:-$_HERE/../2_train/steps}
export ANALYSIS_DIR=${ANALYSIS_DIR:-$S}

say(){ echo "[$(date +%H:%M:%S)] $*"; }

# 코호트 정의 파일 찾기 (COHORT_CONF 주면 그 파일 사용). 검체 수는 여기서 셈
if [ -n "${COHORT_CONF:-}" ]; then
  _CONF=$COHORT_CONF
else
  _CONF=${COHORT_DIR:-$HOME/cohorts}/${COHORT}.conf
fi
[ -f "$_CONF" ] || _CONF=$S/cohorts/${COHORT}.conf
[ -f "$_CONF" ] || { echo "코호트 정의가 없다: ${COHORT}.conf"; exit 1; }
# 1) 검체마다 피처 5종 (run_one_sample 안에서 make_val_config 부름)
# 목록 파일 이름에 코호트와 판 둘 다 넣음 - 코호트 둘 동시 실행 시 충돌 방지
# CRLF 로 저장된 코호트 정의 대비해 CR 제거
_TMP=${TMP_ROOT:-$HOME/tmp}
mkdir -p "$_TMP"
_LIST=$_TMP/samples_${COHORT}_$V.txt
awk -F= '/^SAMPLES=/{sub(/^SAMPLES=/,""); print}' "$_CONF" \
  | tr -d '\r' | tr ' \t' '\n\n' | grep -v '^$' > "$_LIST"
[ -s "$_LIST" ] || { echo "코호트 정의에 SAMPLES 가 비어 있다: $_CONF" >&2; exit 1; }
N=$(wc -l < "$_LIST")
say "검증 시작: 판 $V · 코호트 $COHORT · 검체 $N개 · 동시 $PAR"
say "  코호트 정의  $(cd "$(dirname "$_CONF")" && pwd)/$(basename "$_CONF")"
say "① 피처 만들기"
# NUL 구분(-0) - 검체 이름에 따옴표 있어도 xargs 가 안 멈춤
# xargs 종료값 받아둠 (set -e 는 안 씀, flock 건너뜀에서 죽음)
tr '\n' '\0' < "$_LIST" | xargs -0 -P"$PAR" -I{} bash "$S/run_one_sample.sh" {} "$V"
_rc1=$?
if [ "$_rc1" -ne 0 ]; then
  say "!! ① 에서 실패한 검체가 있습니다 (xargs rc=$_rc1)"
  say "   검체별 로그: ${TMP_ROOT:-$HOME/tmp}/logs/${COHORT}/${V%%_*}/${V}/<검체>.log"
fi

# 2) 검체마다 채점. 필요한 값은 export 말고 위치 인자로 넘김
_APPLY=$S/score_one_sample.py
[ -f "$_APPLY" ] || { echo "score_one_sample.py 가 없습니다: $_APPLY" >&2; exit 1; }
_M=$METH_ROOT/results/panels/${V%%_*}/$V/step04_ML_classifier/models
_nm=$(ls -d "$_M"/*/pooled 2>/dev/null | wc -l)
[ "$_nm" -ge 1 ] || { say "!! 모델이 없습니다 ($V). 2_train 의 09 굳히기를 먼저"; exit 1; }
say "② 채점: pooled 모델 ${_nm}개"
# 채점 로그 남김. 이름에 PID 넣어 드라이버 둘이 같은 파일에 안 씀 (기록용)
_SLOG=${TMP_ROOT:-$HOME/tmp}/logs/${COHORT}/score_${V}_$$.log
mkdir -p "$(dirname "$_SLOG")"
say "  채점 로그  $_SLOG"
tr '\n' '\0' < "$_LIST" | xargs -0 -P"$PAR" -I{} bash -c '
  s=$1; v=$2; apply=$3; root=$4; coh=$5; work=$6; slog=$7; wtmp=$8
  W=$work/${coh}_val/$v/$s
  # 작업 폴더 없으면 1) 이 안 돈 것, 로그에 남기고 넘어감
  [ -d "$W" ] || { echo "[없음] 작업 폴더 $W" >> "$slog"; exit 0; }
  ML=$root/${coh}_val/$s/results/$v/step04_ML_classifier
  # 여기서 건너뛰기 판단 안 함 - score_one_sample.py 가 모델 지문(md5) 보고 직접 판단
  #   셸에서 건너뛰면 모델 다시 굳혀도 옛 점수가 남음
  mkdir -p "$wtmp/lock/val"                           # run_one_sample.sh 와 같은 잠금
  exec 9>"$wtmp/lock/val/${v}__${s}.lock"
  flock -n 9 || { echo "[잠김] 다른 실행이 잡고 있다 $s" >> "$slog"; exit 0; }
  echo "$$ score_one $v $s" >&9
  cd "$W" || { echo "[실패] cd $W" >> "$slog"; exit 1; }
  echo "=== $s $(date +%H:%M:%S)" >> "$slog"
  python -u "$apply" --allmodels >> "$slog" 2>&1 || { echo "!! $s 채점 실패" >> "$slog"; exit 1; }
' _ {} "$V" "$_APPLY" "$METH_ROOT" "$COHORT" "${WORK_ROOT:-$HOME}" "$_SLOG" "$_TMP"
_rc2=$?
[ "$_rc2" -ne 0 ] && say "!! ② 에서 실패한 검체가 있습니다 (xargs rc=$_rc2). 채점 로그 $_SLOG"

# 완료 판정은 로그가 아니라 점수 파일과 잠금으로 함
#   점수 파일 없는 검체는 잠금 잡아 봄 - 잡히면 빠짐, 안 잡히면 다른 실행이 도는 중
# 잠금 폴더 여기서도 만듦 - 없으면 exec 9> 실패가 다른 실행 중으로 잘못 읽힘
mkdir -p "${TMP_ROOT:-$HOME/tmp}/lock/val"
_d=0; _busy=0; _miss=""
while read -r s; do
  # 파일 이름 꼬리 _all_combos, _모든조합 둘 다 완료로 셈
  if [ -s "$METH_ROOT/${COHORT}_val/$s/results/$V/step04_ML_classifier/y_prob_all_${V}_all_combos.csv" ] || [ -s "$METH_ROOT/${COHORT}_val/$s/results/$V/step04_ML_classifier/y_prob_all_${V}_모든조합.csv" ]; then
    _d=$((_d+1)); continue
  fi
  # 잠금 잡았다 바로 놓음. 잡히면 아무도 안 도는 중 (진짜 빠짐)
  if ( exec 9>"${TMP_ROOT:-$HOME/tmp}/lock/val/${V}__${s}.lock"; flock -n 9 ) 2>/dev/null; then
    _miss="$_miss $s"
  else
    _busy=$((_busy+1))
  fi
done < "$_LIST"
say "  채점된 검체 $_d/$N"
# 이번 실행 1), 2) 에서 실패 있으면 중단 (옛 점수 파일 남아 있어도 완료로 안 봄)
if [ "$_rc1" -ne 0 ] || [ "$_rc2" -ne 0 ]; then
  say "!! 이번 실행에서 실패한 검체가 있습니다(① rc=$_rc1 · ② rc=$_rc2). 남은 점수 파일로 완료라 하지 않습니다."
  say "   로그를 보고 고친 뒤 다시 부르십시오. 다음 단계로 넘어가지 마십시오."
  exit 1
fi
say "검증 끝. 검체별 채점: <뿌리>/${COHORT}_val/*/results/$V/step04_ML_classifier/"

# 다 못 끝냈으면 exit 1 (run_validate.sh && make_report.sh 가 빈 표로 안 넘어가게)
# 다른 실행이 잡고 있는 것도 미완료로 셈. 안내 문구만 따로
if [ -n "$_miss" ] || [ "$_busy" -gt 0 ]; then
  [ "$_busy" -gt 0 ] && {
    say "  ${_busy}개는 다른 실행이 「지금」 잡고 있습니다. 그 실행이 끝난 뒤"
    say "  이 명령을 다시 부르면 채워집니다 (고장이 아닙니다)."; }
  [ -n "$_miss" ] && {
    say "!! 빠진 검체:$_miss"
    say "   ① 의 검체별 로그 ${TMP_ROOT:-$HOME/tmp}/logs/${COHORT}/${V%%_*}/${V}/<검체>.log 와"
    say "   ② 의 채점 로그 $_SLOG 을 보십시오."; }
  say "   $_d/$N 만 끝났습니다. 다음 단계로 넘어가지 마십시오."
  exit 1
fi
say "  다음: COHORT=$COHORT bash 4_report/make_report.sh   (pipeline/ 에서 · 표는 위 파일들을 직접 읽는다 · 혈장 코호트는 4_report/steps/plasma_metric.py)"
