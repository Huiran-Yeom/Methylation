#!/bin/bash
# 4단계: 결과 표 정리
#
#   검증 채점 결과 읽어서 표만 출력 (파일 안 씀)
#
#   make_tables.py      판마다 피처 조합 x 비율 AUC 표 (검체 단위, 사람 단위)
#   primary_metric.py   1차 지표: 짝지은 부트스트랩으로 판 사이 차이 구간
#                       음성 고정이라 음성 쪽 불확실성은 구간에 없음 (표에 각주)
#
#   사용법
#     bash make_report.sh                      # 실판 3 + 귀무 6, 5k
#     DEPTH=100k bash make_report.sh           # 깊이를 바꿔서
#     bash make_report.sh j15bl200             # 한 판만
set -u

_HERE=$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)
S=$_HERE/steps

# 공통 머리: config.conf 읽기 (환경변수 우선), BLAS 스레드 고정, 암 쪽 자료 이름 맞추기 (0_setup/load_config.sh)
_LC=""
for _d in "$_HERE/.." "$_HERE/../.."; do
  [ -f "$_d/0_setup/load_config.sh" ] && { _LC="$_d/0_setup/load_config.sh"; break; }
done
: "${_LC:=${METH_CONF_DIR:+$METH_CONF_DIR/load_config.sh}}"
[ -f "${_LC:-}" ] || { echo "멈춤: 0_setup/load_config.sh 를 못 찾았다 — 트리 밖에서 돌리려면 METH_CONF_DIR 를 직접 주라" >&2; exit 1; }
. "$_LC"
# GEN 은 기본값 없음. conf 에 없으면 중단
: "${GEN:?GEN 이 없습니다 — 0_setup/config.conf 의 GEN 을 채우십시오}"
# COHORT 기본값 없음 - 없는 코호트 읽으면 빈 표로 정상 종료돼서
: "${COHORT:?COHORT 를 주십시오 — 3단계에서 쓴 것과 같아야 합니다 (예: COHORT=<your_cohort>)}"
# COHORT 가 VAL_COHORTS 에 있는지 확인 (3단계와 같은 검사)
if [ -n "${VAL_COHORTS:-}" ]; then
  _VC=$(printf %s "$VAL_COHORTS" | tr -d " \t")
  case ",${_VC}," in
    *",${COHORT},"*) : ;;
    *) echo "멈춤: COHORT=$COHORT 는 VAL_COHORTS 에 없습니다 ($VAL_COHORTS)" >&2; exit 1 ;;
  esac
fi
export COHORT
# 혈장 코호트(NEG_SAMPLES 있음)는 이 표 대신 steps/plasma_metric.py 사용
# 코호트 설정 찾는 순서: COHORT_CONF → COHORT_DIR(기본 ~/cohorts) → 저장소 예시
for _c in "${COHORT_CONF:-}" "${COHORT_DIR:-$HOME/cohorts}/$COHORT.conf" "$_HERE/../3_validate/steps/cohorts/$COHORT.conf"; do
  [ -n "$_c" ] && [ -f "$_c" ] || continue
  if grep -qE '^[[:space:]]*NEG_SAMPLES=[[:space:]]*[^[:space:]]' "$_c"; then
    echo "멈춤: $_c 에 NEG_SAMPLES 가 있습니다 — 혈장 코호트는 steps/plasma_metric.py 로 판정합니다" >&2
    exit 1
  fi
  break
done
# DEPTH 는 primary_metric.py 에만 적용 (make_tables.py 는 깊이 전부 돎)
DEPTH=${DEPTH:-5k}

ONE=${1:-}
if [ -n "$ONE" ]; then
  LIST="$ONE"
else
  LIST=""
  _N=${PANEL_N:-200}
  for p in $(echo "${PANELS:-bl$_N,jsd$_N,jsdb$_N}" | tr ',' ' '); do LIST="$LIST j${GEN}${p}"; done
  for p in $(echo "${NULL_PANELS:-}" | tr ',' ' '); do LIST="$LIST j${GEN}${p}"; done
fi

echo "=== 판별 표 (코호트 $COHORT · 깊이 $DEPTH) ==="
_bad=0
for p in $LIST; do
  echo
  echo "##### $p"
  # 실패한 판 기록 (끝에서 종료값으로 알림)
  ( cd "$S" && python -u make_tables.py "$p" ) || { echo "!! 표 실패 $p" >&2; _bad=1; }
done

echo
echo "=== 1차 지표: 짝지은 부트스트랩 (깊이 $DEPTH) ==="
( cd "$S" && python -u primary_metric.py --depth "$DEPTH" ) || _bad=1

# 하나라도 실패하면 rc 1
if [ "$_bad" -ne 0 ]; then
  echo "!! 실패한 표가 있습니다. 위 !! 줄을 보십시오" >&2
  exit 1
fi
