#!/bin/bash
# ── 3단계 · 외부검증 ──────────────────────────────────────────────────
#
#   학습에 쓰지 않은 정상을 배경으로, 외부 코호트의 암 검체를 학습과 「똑같은」
#   비율·깊이·믹싱 방식으로 섞어 학습된 모델로 채점한다.
#   바뀌는 것은 입력 BAM 과 출력 경로뿐이다. 계산은 학습과 같아야 한다.
#
#   ★ 이 단계는 학습 스크립트를 검체 폴더로 「복사해」 돌린다.
#     그래서 2_train/steps 의 파일 이름이 config.conf 의 STEP_MAP 과 엮여 있다.
#     이름을 바꾸면 STEP_MAP 도 같이 고친다.
#
#   흐름
#     검체마다   make_val_config -> run_one_sample (피처 다섯) -> score_one_sample (채점)
#     끝난 뒤    4_report/make_report.sh 가 검체별 채점 파일을 직접 읽어 표를 만든다
#
#   사용
#     bash run_validate.sh <판이름> [동시실행수]
#     COHORT=<your_cohort> bash run_validate.sh j15bl200_5k_uni 8
set -u

_HERE=$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)
S=$_HERE/steps

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
: "${COHORT:?COHORT 를 주십시오 (예: COHORT=<your_cohort>)}"

# 2026-09-29: flock 이 없으면 아래 회계가 「잠금을 못 잡았다」 로 읽어 전부
#   「다른 실행이 쥐고 있다」가 된다. 판정이 거꾸로 선다. 먼저 확인한다.
#   (껍데기 시험에서 실제로 그렇게 나왔다. Git Bash 에는 flock 이 없다.)
command -v flock >/dev/null 2>&1 || {
  echo "멈춤: flock 이 없습니다 (util-linux). 이 파이프라인은 검체별 잠금에 씁니다." >&2
  echo "      없으면 「빠짐」과 「남이 쥐고 있음」을 구분하지 못합니다." >&2
  exit 1; }

# 2026-09-29: 예전에는 이 관문이 aggregate_scores.py 에 있었다. 그 파일을 이 저장소에서
#   빼면서 관문도 같이 사라졌으므로, 코호트를 「고르는 자리」 인 여기로 옮긴다.
#   목적: 옛 코호트나 오타를 넣으면 조용히 빈 결과를 내는 대신 멈춘다.
#   VAL_COHORTS 가 비어 있으면 아무 코호트나 허용한다 (내려받은 그대로 쓸 수 있게).
#   구분자 뒤 공백을 허용한다. `a, b` 로 적으면 b 를 「거절」 했다(2026-09-29).
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
# 검증이 학습 스크립트를 복사해 오는 자리와 make_val_config 위치
export TRAIN_SRC=${TRAIN_SRC:-$_HERE/../2_train/steps}
export ANALYSIS_DIR=${ANALYSIS_DIR:-$S}

say(){ echo "[$(date +%H:%M:%S)] $*"; }

# 코호트 정의에서 검체 수를 센다. 박지 않는다
# 2026-09-29: COHORT_CONF 로 못 박을 수 있게 한다. 읽은 자리를 절대경로로 찍는다.
if [ -n "${COHORT_CONF:-}" ]; then
  _CONF=$COHORT_CONF
else
  _CONF=${COHORT_DIR:-$HOME/cohorts}/${COHORT}.conf
fi
[ -f "$_CONF" ] || _CONF=$S/cohorts/${COHORT}.conf
[ -f "$_CONF" ] || { echo "코호트 정의가 없다: ${COHORT}.conf"; exit 1; }
# ① 검체마다 피처 다섯 (run_one_sample 안에서 make_val_config 를 부른다)
#   목록 이름에 코호트와 판을 「둘 다」 넣는다. 같은 판을 코호트 둘로 동시에 돌리면
#   판만으로는 서로의 목록을 밟는다(2026-09-29).
#   `tr -d '\r'` : 코호트 정의가 CRLF 로 저장돼 있으면 마지막 검체 이름에 CR 이
#   붙어 경로가 어긋나고, 그런데도 wc -w 는 그것을 하나로 세어 「1개 빠짐」으로만
#   보인다. 여기서 한 번 걷어낸다. [CRLF 가 셸을 죽인다]
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
# -0 : 검체 이름에 따옴표가 들어가면 xargs -I{} 는 「전체 묶음」 을 unmatched quote 로
#   중단한다. 뒤의 검체가 조용히 안 돈다. NUL 구분이면 그 부류가 없다.
# 2026-09-29: xargs 의 종료값을 버리지 않는다. run_one_sample.sh 는 원본 폴더를 못
#   찾거나 단계 파일이 없으면 exit 1 하는데(둘 다 AUC 0.500 길목이다), 그 신호가
#   여기서 사라지면 「검체별 로그를 열어 보지 않는 한」 아무 표시가 없다.
#   이 스크립트는 set -e 를 쓰지 않으므로 종료값만 받아 두면 된다 (set -e 를 켜면
#   flock 건너뜀 같은 「정상적인 0 아님」 에서 죽는다).
tr '\n' '\0' < "$_LIST" | xargs -0 -P"$PAR" -I{} bash "$S/run_one_sample.sh" {} "$V"
_rc1=$?
if [ "$_rc1" -ne 0 ]; then
  say "!! ① 에서 실패한 검체가 있습니다 (xargs rc=$_rc1)"
  say "   검체별 로그: ${TMP_ROOT:-$HOME/tmp}/logs/${COHORT}/${V%%_*}/${V}/<검체>.log"
fi

# ② 검체마다 채점
#   2026-09-29: 예전에는 별도 병렬 런처(score_parallel.sh)가 임시 스크립트를 써서
#     돌렸다. 그 방식에 결함이 둘 있었다. (가) 임시 스크립트가 단일인용 heredoc
#     안에서 $_APPLY 를 쓰는데 export 하지 않아 검체 전부가 조용히 실패할 수 있었고,
#     (나) 검체 목록을 valj_samples.txt 로 읽는데 여기서는 valj_samples_$V.txt 를
#     쓰고 있어 아예 다른 파일을 봤다.
#   그래서 런처를 없애고 여기서 바로 돈다. 필요한 값은 「위치 인자로」 넘긴다 —
#     export 에 기대지 않으면 (가) 부류가 구조적으로 생기지 않는다.
_APPLY=$S/score_one_sample.py
[ -f "$_APPLY" ] || { echo "score_one_sample.py 가 없습니다: $_APPLY" >&2; exit 1; }
_M=$METH_ROOT/results/panels/${V%%_*}/$V/step04_ML_classifier/models
_nm=$(ls -d "$_M"/*/pooled 2>/dev/null | wc -l)
[ "$_nm" -ge 1 ] || { say "!! 모델이 없습니다 ($V) — 2_train 의 09 굳히기를 먼저"; exit 1; }
say "② 채점: pooled 모델 ${_nm}개"
# 채점 로그를 버리지 않는다. 예전에는 > /dev/null 2>&1 이라 결측 특징 %·비율 건너뜀·
#   트레이스백이 전부 사라지고 「파일이 있나」만 남았다(2026-09-29).
# 2026-09-29: 이름에 PID 를 넣는다. 코호트·판으로만 갈랐더니 드라이버 둘이 같은
#   파일에 썼고(run_one_sample.sh:73 이 같은 이유로 코호트·판을 다 넣는다),
#   지난 줄이 섞여 판정을 흐렸다. 이제 판정은 잠금 실측으로 하므로 로그는 기록용이다.
_SLOG=${TMP_ROOT:-$HOME/tmp}/logs/${COHORT}/score_${V}_$$.log
mkdir -p "$(dirname "$_SLOG")"
say "  채점 로그  $_SLOG"
tr '\n' '\0' < "$_LIST" | xargs -0 -P"$PAR" -I{} bash -c '
  s=$1; v=$2; apply=$3; root=$4; coh=$5; work=$6; slog=$7; wtmp=$8
  W=$work/${coh}_val/$v/$s
  # 작업 폴더가 없으면 ① 이 안 돈 것이다. 조용히 넘기면 「전부 건너뜀」 이
  #   0/N 으로만 보인다. 무엇이 없었는지 남긴다.
  [ -d "$W" ] || { echo "[없음] 작업 폴더 $W" >> "$slog"; exit 0; }
  ML=$root/${coh}_val/$s/results/$v/step04_ML_classifier
  # 2026-09-29: 여기서 「model_md5 가 있으면 건너뛴다」를 「하지 않는다」.
  #   score_one_sample.py 의 _already_done() 은 표의 지문을 「현재」 BEST.joblib 의
  #   md5 와 맞춰 보고, 모델을 다시 굳히면 스스로 다시 채점한다. 셸에서 헤더만
  #   보고 건너뛰면 그 자기치유가 영구히 죽어 옛 모델 점수가 표에 남는다.
  mkdir -p "$wtmp/lock/val"                           # run_one_sample.sh 와 같은 잠금
  exec 9>"$wtmp/lock/val/${v}__${s}.lock"
  flock -n 9 || { echo "[잠김] 다른 실행이 잡고 있다 $s" >> "$slog"; exit 0; }
  echo "$$ score_one $v $s" >&9
  cd "$W" || { echo "[실패] cd $W" >> "$slog"; exit 1; }
  echo "=== $s $(date +%H:%M:%S)" >> "$slog"
  python -u "$apply" --allmodels >> "$slog" 2>&1 || echo "!! $s 채점 실패" >> "$slog"
' _ {} "$V" "$_APPLY" "$METH_ROOT" "$COHORT" "${WORK_ROOT:-$HOME}" "$_SLOG" "$_TMP"

# 2026-09-29: 판정 근거를 「로그 세기」 에서 「잠금 실측」 으로 바꾼다.
#   로그로 세면 셋이 어긋난다. (가) 「잠겨 건너뜀」은 ② 에만 적히는데 이어달리기
#   충돌은 ① 에서 나고, (나) 로그 파일을 코호트·판으로만 갈라 드라이버 둘이 같은
#   파일을 쓰며, (다) 「파일 있음」과 「이번에 건너뜀」이 겹쳐 두 번 세어진다.
#   그래서 남은 검체마다 「그 잠금을 지금 누가 쥐고 있나」 를 직접 본다. 겹치지 않는다.
# 2026-09-29 ②: 잠금 폴더를 「여기서도」 만든다. 없으면 `exec 9>` 가 flock 에 닿기도
#   전에 실패하고, 그 실패가 「다른 실행이 쥐고 있다」로 읽힌다. 즉 ① 이 전 검체에서
#   죽은 상황이 「실패 아님 · exit 0」 이 된다. 판정이 정확히 거꾸로 선다.
mkdir -p "${TMP_ROOT:-$HOME/tmp}/lock/val"
_d=0; _busy=0; _miss=""
while read -r s; do
  if [ -s "$METH_ROOT/${COHORT}_val/$s/results/$V/step04_ML_classifier/y_prob_all_${V}_all_combos.csv" ]; then
    _d=$((_d+1)); continue
  fi
  # 잠금을 잡아 보고 곧 놓는다. 잡히면 아무도 안 하고 있다(진짜 빠짐).
  if ( exec 9>"${TMP_ROOT:-$HOME/tmp}/lock/val/${V}__${s}.lock"; flock -n 9 ) 2>/dev/null; then
    _miss="$_miss $s"
  else
    _busy=$((_busy+1))
  fi
done < "$_LIST"
say "  채점된 검체 $_d/$N"
say "검증 끝 — 검체별 채점: <뿌리>/${COHORT}_val/*/results/$V/step04_ML_classifier/"

# 다 못 했으면 「0 이 아닌 값」 으로 끝낸다. 예전에는 마지막 명령이 say 라 0/N 이어도
#   exit 0 이었고, `run_validate.sh && make_report.sh` 가 빈 표로 넘어갔다.
# 2026-09-29 ②: 「다른 실행이 쥐고 있다」도 「이 실행은 못 끝냈다」 이다. 뒤 단계로
#   넘어가면 안 된다. 표는 빠진 것을 메우지 않는다. 다만 원인이 다르므로 말은 나눈다.
if [ -n "$_miss" ] || [ "$_busy" -gt 0 ]; then
  [ "$_busy" -gt 0 ] && {
    say "  ${_busy}개는 다른 실행이 「지금」 잡고 있습니다 — 그 실행이 끝난 뒤"
    say "  이 명령을 다시 부르면 채워집니다 (고장이 아닙니다)."; }
  [ -n "$_miss" ] && {
    say "!! 빠진 검체:$_miss"
    say "   ① 의 검체별 로그 ${TMP_ROOT:-$HOME/tmp}/logs/${COHORT}/${V%%_*}/${V}/<검체>.log 와"
    say "   ② 의 채점 로그 $_SLOG 을 보십시오."; }
  say "   $_d/$N 만 끝났습니다. 다음 단계로 넘어가지 마십시오."
  exit 1
fi
say "  다음: bash ../4_report/make_report.sh   (표는 위 파일들을 직접 읽는다)"
