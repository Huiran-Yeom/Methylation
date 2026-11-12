#!/bin/bash
# 1단계: DMR 패널 선정
#
#   후보 블록 → 창별 무늬 분포로 JSD → 패널 PANEL_N 블록 (기본 200)
#   비교용 무작위 패널 6개도 여기서 만듦
#
#   설정은 ../0_setup/config.conf, 단계 코드는 steps/ (번호 = 실행 순서)
#
#   사용법
#     bash run_dmr.sh              # 전 단계
#     bash run_dmr.sh 04           # 한 단계만
#     JSD_LIMIT=500000 bash run_dmr.sh   # 리드 줄여서 도는지만 확인 (값은 다름)
#
#   오래 걸림: 01 (cov 전체), 03 (BAM 전체), 07 (부트스트랩 100회)
set -u

_HERE=$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)
S=$_HERE/steps
ONE=${1:-}

# 공통 머리: config.conf 읽기 (환경변수 우선), BLAS 스레드 고정, 암 쪽 자료 이름 맞추기 (0_setup/load_config.sh)
_LC=""
for _d in "$_HERE/.." "$_HERE/../.."; do
  [ -f "$_d/0_setup/load_config.sh" ] && { _LC="$_d/0_setup/load_config.sh"; break; }
done
: "${_LC:=${METH_CONF_DIR:+$METH_CONF_DIR/load_config.sh}}"
[ -f "${_LC:-}" ] || { echo "멈춤: 0_setup/load_config.sh 를 못 찾았다. 트리 밖에서 돌리려면 METH_CONF_DIR 를 직접 주라" >&2; exit 1; }
. "$_LC"
: "${METH_ROOT:=/ssd_data/Methylation}"
# GEN 은 기본값 없음. conf 에 없으면 중단
: "${GEN:?GEN 이 없습니다. 0_setup/config.conf 의 GEN 을 채우십시오}"
# 이 단계 경로. 폴더 구조 다르면 여기만 수정
D=$METH_ROOT/results/dmr                              # 1단계 산출 폴더
ND=${NORMAL_SET:-$METH_ROOT/sample_data/_input/normal_pub15}   # 정상 cov (*.cov.gz)

say(){ echo "[$(date +%H:%M:%S)] $*"; }
# 산출 폴더 비어있지 않으면 중단 (덮어쓰려면 OVERWRITE=<폴더명>)
_check() {   # _check <폴더>
  local d=$1
  [ -d "$d" ] || return 0
  [ -z "$(ls -A "$d" 2>/dev/null)" ] && return 0
  if [ "${OVERWRITE:-}" = "$(basename "$d")" ]; then
    echo "[덮어쓴다] $d  (OVERWRITE=$(basename "$d") 를 받았다)" >&2
    return 0
  fi
  echo "" >&2
  echo "멈춤: 산출 폴더가 비어 있지 않습니다. $d" >&2
  echo "  덮어쓰면 기존 결과가 사라집니다. 셋 중 하나를 고르십시오:" >&2
  echo "    1) 다른 뿌리에 쓴다   METH_ROOT=\$HOME/tmp/test bash $(basename "$0")" >&2
  echo "    2) 새 세대로 간다     GEN=16 bash $(basename "$0")" >&2
  echo "    3) 정말 덮는다        OVERWRITE=$(basename "$d") bash $(basename "$0")" >&2
  echo "" >&2
  exit 1
}

# 단계별로 새로 만드는 산출 폴더 (한 단계만 돌릴 때도 그 폴더만 검사)
_outs() {   # _outs <단계>
  # 앞 단계 폴더에 추가만 하는 단계(02, 04b, 04c, 09, 10, 10c)는 빈 값
  # 전체 실행이면 앞 단계에서 이미 검사됨
  case "$1" in
    01)  echo "$D/j_candidates$GEN" ;;
    02)  : ;;                       # j_candidates 에 추가
    03)  echo "$D/j_jsdcount$GEN" ;;
    04)  echo "$D/j_jsd$GEN" ;;
    04b|04c) : ;;                   # j_jsd 에 추가
    05)  echo "$D/j_llrref$GEN" ;;
    06)  echo "$D/j_panel$GEN" ;;
    06b) echo "$D/j_panel${GEN}_supplement" ;;
    07)  echo "$D/j_panel_dmr$GEN/bl${PANEL_N}_d10_cellline" ;;
    08)  echo "$D/j_panel_dmr${GEN}_panels" ;;
    09|10) : ;;                     # j_panel_dmr<GEN>_panels 에 추가
    10c) : ;;                      # j_panel<GEN> 에 추가
    11)  echo "$D/j_pool${GEN}_cellline"     # 풀 3개 다 검사
         echo "$D/j_pool${GEN}_rand"
         echo "$D/j_pool${GEN}_randb" ;;
    *)   : ;;
  esac
}

step(){ # step <번호> <제목> <명령...>
  local no=$1 title=$2; shift 2
  if [ -n "$ONE" ] && [ "$ONE" != "$no" ]; then return 0; fi
  say "── $no  $title"
  "$@" || { say "!! 실패 $no"; exit 1; }
}

# 실행 전에 이번에 돌 단계들의 산출 폴더 한 번에 검사
# 단계 목록은 이 파일의 step 줄에서 읽음
_precheck() {
  local _ids _id _o _seen=" "
  if [ -n "$ONE" ]; then
    _ids=$ONE
  else
    # 목록 못 읽으면(파이프로 넘긴 경우 등) 검사 없이 지나가지 않게 중단
    _ids=$(grep -oE "^step +[0-9a-z]+" "$0" 2>/dev/null | awk "{print \$2}")
    if [ -z "$_ids" ]; then
      echo "멈춤: 단계 목록을 못 읽어 중단했습니다 ($0). 이대로 돌면 덮어쓰기 관문이 꺼집니다." >&2
      echo "      파이프로 넘기지 말고 파일로 저장한 뒤 bash $(basename "$0") 로 실행하십시오." >&2
      exit 1
    fi
  fi
  # 없는 단계 번호면 아무것도 안 하고 끝나버리니까 먼저 중단
  if [ -n "$ONE" ]; then
    case " $(grep -oE "^step +[0-9a-z]+" "$0" 2>/dev/null | awk "{print \$2}" | tr "
" " ") " in
      *" $ONE "*) : ;;
      *) echo "멈춤: 그런 단계가 없습니다. $ONE" >&2
         echo "  있는 단계: $(grep -oE "^step +[0-9a-z]+" "$0" | awk "{print \$2}" | tr "
" " ")" >&2
         exit 1 ;;
    esac
  fi
  for _id in $_ids; do
    while IFS= read -r _o; do
      [ -n "$_o" ] || continue
      case "$_seen" in *"|$_o|"*) continue ;; esac
      _seen="$_seen|$_o|"
      _check "$_o"
    done <<EOF
$(_outs "$_id")
EOF
  done
}

# 패널 크기. PANELS, NULL_PANELS 이름 속 크기도 같아야 함 (06 이후 단계만 봄, 01~05 는 크기와 무관)
: "${PANEL_N:=200}"
case "$PANEL_N" in ''|*[!0-9]*|0) echo "멈춤: PANEL_N 은 1 이상의 정수여야 합니다: $PANEL_N" >&2; exit 1 ;; esac
_TC=${TOPCAP:-1000}; [ "$_TC" = N ] && _TC=$PANEL_N      # 07 --topcap. N 이면 PANEL_N 과 같게
_names_ok() {
  local x bad=""
  for x in ${PANELS//,/ }; do
    [[ $x =~ ^(bl|jsd|jsdb)([0-9]+)$ ]] && [ "${BASH_REMATCH[2]}" = "$PANEL_N" ] || bad="$bad $x"
  done
  for x in ${NULL_PANELS//,/ }; do
    [[ $x =~ ^(randb?)([1-3])([0-9]+)$ ]] && [ "${BASH_REMATCH[3]}" = "$PANEL_N" ] || bad="$bad $x"
  done
  [ -z "$bad" ] && return 0
  echo "멈춤: PANELS · NULL_PANELS 이름이 PANEL_N($PANEL_N) 과 안 맞습니다:$bad" >&2
  echo "  N 을 바꾸면 config.conf 의 PANELS, NULL_PANELS 이름 숫자도 같이 바꾸십시오 (예 bl300, rand1300)" >&2
  exit 1
}

say "DMR 시작: METH_ROOT=$METH_ROOT · GEN=$GEN"
case "${ONE:-all}" in 0[1-5]|04b|04c) ;; *) _names_ok ;; esac
_precheck

step 01 "후보 블록"        python -u "$S/01_candidates.py" --nd "$ND" --out "$D/j_candidates$GEN"
step 02 "블록별 CpG"       python -u "$S/02_blockcpg.py"  --nd "$ND" --cand "$D/j_candidates$GEN"
# JSD_LIMIT: 리드 수 제한 (도는지 확인용, 본 실행에선 안 씀)
_LIM=""
[ -n "${JSD_LIMIT:-}" ] && _LIM="--limit $JSD_LIMIT"
step 03 "창별 무늬수"      python -u "$S/03_jsdcount.py"  --nd "$ND" --cand "$D/j_candidates$GEN" \
                                --out "$D/j_jsdcount$GEN" $_LIM
step 04 "JSD"              python -u "$S/04_jsd.py"       --dir "$D/j_jsdcount$GEN" \
                                --cand "$D/j_candidates$GEN" --out "$D/j_jsd$GEN" --k "${PATTERN_K:-3}" --minr "${MIN_READS_PER_WINDOW:-6}"
# 04b(귀무 분포), 04c(부트스트랩)는 06b 의 jsdb<N> 선정에 필요
step 04b "JSD 귀무 분포"    python -u "$S/04b_jsdnull.py"  --dir "$D/j_jsdcount$GEN" \
                                --cand "$D/j_candidates$GEN" --out "$D/j_jsd$GEN" --k "${PATTERN_K:-3}" --minr "${MIN_READS_PER_WINDOW:-6}"
step 04c "JSD 부트스트랩"   python -u "$S/04c_jsdboot.py"  --dir "$D/j_jsdcount$GEN" \
                                --cand "$D/j_candidates$GEN" --out "$D/j_jsd$GEN" --k "${PATTERN_K:-3}" --minr "${MIN_READS_PER_WINDOW:-6}"
step 05 "LLR 참조표"       python -u "$S/05_llrref.py"    --dir "$D/j_jsdcount$GEN" \
                                --cand "$D/j_candidates$GEN" --out "$D/j_llrref$GEN" --k "${PATTERN_K:-3}" --minr "${MIN_READS_PER_WINDOW:-6}"

_JF=$(ls "$D/j_jsd$GEN"/jsd_K${PATTERN_K:-3}_minr${MIN_READS_PER_WINDOW:-6}.parquet 2>/dev/null | head -1)
step 06 "JSD 상위 $PANEL_N 선정" python -u "$S/06_select.py" --cand "$D/j_candidates$GEN" \
                                --jsd "$_JF" --out "$D/j_panel$GEN" --n "$PANEL_N"
# 06b 는 PANELS 에 jsdb<N> 있을 때만. N 블록 못 채우는 자료면 PANELS 에서 빼기
_06b() {
  case ",${PANELS:-bl$PANEL_N,jsd$PANEL_N,jsdb$PANEL_N}," in
    *,jsdb$PANEL_N,*) python -u "$S/06b_select_supplement.py" ;;
    *) echo "  건너뜀. PANELS(${PANELS}) 에 jsdb$PANEL_N 이 없다" ;;
  esac
}
step 06b "보완판 JSD $PANEL_N (jsdb$PANEL_N)" _06b
step 07 "Baseline 순위 $PANEL_N" python -u "$S/07_panel.py"  --cand "$D/j_candidates$GEN" \
                                --out "$D/j_panel_dmr$GEN/bl${PANEL_N}_d10_cellline" --topn "$PANEL_N" \
                                --topcap "$_TC"

step 08 "검증용 형식으로 내보내기" python -u "$S/08_export_panels.py"
step 09 "귀무A 3판 (씨앗 20260914+i)"      python -u "$S/09_null_a.py"
step 10 "귀무B 3판 (씨앗 20260914+100+i)" python -u "$S/10_null_b.py"

# 11: 2단계 섞기가 읽는 풀 BAM 생성
# 풀은 자기 패널 영역 리드만 들어있어서 cellline / rand / randb 따로 만듦
# 입력 목록은 10c 가 만듦. USE_NORMAL, USE_GBM 주면 그걸 사용
: "${USE_NORMAL:=$D/j_panel$GEN/use${GEN}_normal.txt}"
: "${USE_GBM:=$D/j_panel$GEN/use${GEN}_gbm.txt}"
export USE_NORMAL USE_GBM

_pool_all() {   # 풀 3개 순서대로, 하나라도 실패하면 중단
  local w bad=""
  for w in cellline rand randb; do
    say "   풀 $w"
    bash "$S/run_pool.sh" "$w" || bad="$bad $w"
  done
  [ -z "$bad" ] && return 0
  echo "멈춤: 풀을 못 만들었습니다:$bad" >&2
  echo "  cellline 은 bl·jsd·jsdb 실판이, rand 는 rand1-3 이," >&2
  echo "  randb 는 randb1-3 이 씁니다. 빠진 풀의 판은 학습이 $PANEL_N 칸 중" >&2
  echo "  몇 칸으로 조용히 돌거나 「풀 부족」 으로 멈춥니다." >&2
  return 1
}
step 10c "풀 입력 만들기 (합집합 3 · 검체목록 2)" python -u "$S/10c_panel_union.py"

step 11 "패널 영역 풀 BAM (cellline·rand·randb)" _pool_all

say "DMR 끝: $D/j_panel_dmr${GEN}_panels/<패널>/DMR_confirmed_<패널>.csv"
