#!/bin/bash
# ── 4단계 · 결과 정리 ─────────────────────────────────────────────────
#
#   검증 채점 결과를 읽어 논문에 실을 표를 만든다. 읽기만 한다 — 아무것도 안 쓴다.
#
#   make_tables.py    판마다 «피처 조합 × 비율» AUC 표. 검체 단위와 사람 단위를 같이 낸다.
#   primary_metric.py   사전등록 1차 지표 — 짝지은 부트스트랩으로 판 사이 차이의 구간.
#                음성을 고정하므로 구간에 «음성 축 불확실성이 0 으로» 들어간다.
#                각주로 표에 같이 박는다.
#
#   사용
#     bash make_report.sh                      # 실판 3 + 귀무 6, 5k
#     DEPTH=100k bash make_report.sh           # 깊이를 바꿔서
#     bash make_report.sh j15bl200             # 한 판만
set -u

_HERE=$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)
S=$_HERE/steps

# 2026-09-28: 보관 대상을 손으로 들고 있으면 목록에 없는 키는 conf 가 환경변수를
#   덮는다 — meth_config.py 의 「환경변수가 이긴다」와 반대가 된다. conf 에 적힌 키를
#   «전부» 보관한다. 목록을 유지할 필요가 없어진다.
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
# 2026-09-28: 셸이 기본값을 채우면 conf 를 비웠을 때 관문이 못 잔다.
#   정본은 config.conf 와 meth_config.py 둘이다. 셸은 값을 만들지 않는다.
: "${GEN:?GEN 이 없습니다 — 0_setup/config.conf 의 GEN 을 채우십시오}"
# 2026-09-29: 코호트 기본값을 우리 내부 이름으로 박지 않는다. 남이 받아서
#   COHORT 를 안 주면 조용히 «없는 코호트» 를 읽고 표 아홉 개가 다 비는데
#   exit 0 으로 끝났다. 주지 않으면 멈춘다.
: "${COHORT:?COHORT 를 주십시오 — 3단계에서 쓴 것과 같아야 합니다 (예: COHORT=<your_cohort>)}"
# 3단계와 같은 관문을 여기에도 둔다. 코호트를 «고르는 자리» 가 둘이면 관문도 둘이어야
#   한다 — 3단계만 막으면 오타난 코호트로 표를 뽑는 길이 열린 채 남는다.
if [ -n "${VAL_COHORTS:-}" ]; then
  _VC=$(printf %s "$VAL_COHORTS" | tr -d " \t")
  case ",${_VC}," in
    *",${COHORT},"*) : ;;
    *) echo "멈춤: COHORT=$COHORT 는 VAL_COHORTS 에 없습니다 ($VAL_COHORTS)" >&2; exit 1 ;;
  esac
fi
export COHORT
# 2026-09-29: DEPTH 는 primary_metric.py --depth 에만 간다. make_tables.py 는
#   깊이를 스스로 돌기 때문에 이 값에 영향받지 않는다 — 머리글이 「표가 바뀐다」로
#   읽히지 않게 적어 둔다.
DEPTH=${DEPTH:-5k}

ONE=${1:-}
if [ -n "$ONE" ]; then
  LIST="$ONE"
else
  LIST=""
  for p in $(echo "${PANELS:-bl200,jsd200,jsdb200}" | tr ',' ' '); do LIST="$LIST j${GEN}${p}"; done
  for p in $(echo "${NULL_PANELS:-}" | tr ',' ' '); do LIST="$LIST j${GEN}${p}"; done
fi

echo "=== 판별 표 (코호트 $COHORT · 깊이 $DEPTH) ==="
_bad=0
for p in $LIST; do
  echo
  echo "##### $p"
  # 2026-09-29: 종료값을 버리면 표 아홉 중 하나가 죽어도 다른 여덟에 묻힌다.
  ( cd "$S" && python -u make_tables.py "$p" ) || { echo "!! 표 실패 $p" >&2; _bad=1; }
done

echo
echo "=== 1차 지표 — 짝지은 부트스트랩 (깊이 $DEPTH) ==="
( cd "$S" && python -u primary_metric.py --depth "$DEPTH" ) || _bad=1

# 2026-09-29: 하나라도 실패했으면 0 이 아닌 값으로 끝낸다. 표 아홉 중 하나가
#   죽어도 나머지 출력에 묻혀 «성공» 으로 보이던 것을 막는다.
if [ "$_bad" -ne 0 ]; then
  echo "!! 실패한 표가 있습니다 — 위 !! 줄을 보십시오" >&2
  exit 1
fi
