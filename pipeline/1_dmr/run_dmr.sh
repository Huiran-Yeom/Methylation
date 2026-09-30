#!/bin/bash
# ── 1단계 · DMR 패널 선정 ─────────────────────────────────────────────
#
#   후보 블록을 고르고, 창별 무늬 분포로 JSD 를 재서 패널 200칸을 확정한다.
#   귀무 패널 6판(무작위 추출)도 여기서 만든다. 되돌림 조건 ㉠ 의 실측 잣대다.
#
#   설정은 ../0_setup/config.conf 한 곳에서만 고친다.
#   각 단계의 실제 코드는 steps/ 에 있다. 번호가 곧 실행 순서다.
#
#   사용
#     bash run_dmr.sh              # 전 단계
#     bash run_dmr.sh 04           # 그 단계만
#     JSD_LIMIT=500000 bash run_dmr.sh   # 리드를 제한해 「도는지」 만 확인 (값은 실제와 다르다)
#
#   오래 걸리는 것: 01(cov 전수 훑기) · 03(BAM 전수 훑기) · 07(부트스트랩 100회)
set -u

_HERE=$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)
S=$_HERE/steps
ONE=${1:-}

# 환경변수로 주는 값이 이긴다 (0_setup/meth_config.py 의 값() 과 같은 순서).
# 보관 대상을 손으로 들고 있으면 목록에 없는 키는 conf 가 환경변수를
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
: "${METH_ROOT:=/ssd_data/Methylation}"
# 셸이 기본값을 채우면 conf 를 비웠을 때 관문이 못 잔다.
#   정본은 config.conf 와 meth_config.py 둘이다. 셸은 값을 만들지 않는다.
: "${GEN:?GEN 이 없습니다 — 0_setup/config.conf 의 GEN 을 채우십시오}"
# ╔══ 이 단계가 쓰는 경로 ═══════════════════════════════════════════╗
#   배치가 다르면 「여기만」 고친다. config.conf 에는 뿌리 둘과 실험 설정만 둔다.
#   정상 cov 이름이 여기(pub15)와 스텝 기본값(pub20)에서 달랐다 —
#     드라이버가 --nd 로 덮어써서 가려져 있었다. 이제 스텝도 pub15 를 쓴다.
D=$METH_ROOT/results/dmr                              # 1단계 산출물이 쌓이는 곳
ND=${NORMAL_SET:-$METH_ROOT/sample_data/_input/normal_pub15}   # 정상 cov (*.cov.gz)
# ╚══════════════════════════════════════════════════════════════════╝

say(){ echo "[$(date +%H:%M:%S)] $*"; }
# ── [데이터가 막는다] ────────────────────────────────────
#   어제 격리 없이 이 드라이버를 돌려 실제 산출물을 덮을 뻔했다(30초 만에 멈춰 피해 0).
#   관문들은 이걸 못 막는다. 드라이버를 「정상 경로」 로 부른 것이니까.
#   플래그(--confirm)는 일주일이면 손이 외워서 관문이 아니게 되고, 「격리는 부르는
#   쪽 책임」은 방심한 사람이 부를 때 정확히 실패한다.
#   그래서 판정 근거를 데이터에 둔다: 대상 폴더가 비어 있지 않으면 멈춘다.
#   덮어야 하면 폴더 「이름을 다시 쳐서」 OVERWRITE 로 준다 — 관성으로 넘길 수 없다.
_check() {   # _check <폴더>
  local d=$1
  [ -d "$d" ] || return 0
  [ -z "$(ls -A "$d" 2>/dev/null)" ] && return 0
  if [ "${OVERWRITE:-}" = "$(basename "$d")" ]; then
    echo "[덮어쓴다] $d  (OVERWRITE=$(basename "$d") 를 받았다)" >&2
    return 0
  fi
  echo "" >&2
  echo "멈춤: 산출 폴더가 비어 있지 않습니다 — $d" >&2
  echo "  덮어쓰면 기존 결과가 사라집니다. 셋 중 하나를 고르십시오:" >&2
  echo "    1) 다른 뿌리에 쓴다   METH_ROOT=\$HOME/tmp/test bash $(basename "$0")" >&2
  echo "    2) 새 세대로 간다     GEN=16 bash $(basename "$0")" >&2
  echo "    3) 정말 덮는다        OVERWRITE=$(basename "$d") bash $(basename "$0")" >&2
  echo "" >&2
  exit 1
}

# 단계마다 「그 단계가 쓰는 폴더」 를 적어 둔다.
#   앞서 한 번은 후보 폴더 하나만 봤고(완성된 패널·풀이 조용히 덮였다), 그 다음엔
#   단일 단계 실행일 때 관문을 통째로 껐다. 그건 더 나쁘다. `bash run_dmr.sh 01` 이
#   보호 없이 j_candidates 를 덮게 되고, 그게 실수로 덮는 가장 흔한 길이다.
#   단계가 쓰는 곳만 보면 우회도 오탐도 없다.
_outs() {   # _outs <단계번호> -> 그 단계가 「처음 만드는」 폴더
  #   앞 단계가 이미 쓴 폴더에 「더하는」 단계(02·04b·04c·09·10)는 빈 값을 돌려준다.
  #   그러지 않으면 `bash run_dmr.sh 04b` 가 04 의 산출을 보고 늘 멈춘다 —
  #   04b 는 그 폴더가 「차 있어야」 도는 단계다. 전체 실행에서는 앞 단계(01·04·08)가
  #   같은 폴더를 이미 검사하므로 보호에 구멍이 생기지 않는다.
  case "$1" in
    01)  echo "$D/j_candidates$GEN" ;;
    02)  [ -n "$ONE" ] && echo "$D/j_candidates$GEN" || : ;;   # 단독 실행이면 검사
    03)  echo "$D/j_jsdcount$GEN" ;;
    04)  echo "$D/j_jsd$GEN" ;;
    04b|04c) [ -n "$ONE" ] && echo "$D/j_jsd$GEN" || : ;;
    05)  echo "$D/j_llrref$GEN" ;;
    06)  echo "$D/j_panel$GEN" ;;
    06b) echo "$D/j_panel${GEN}_supplement" ;;
    07)  echo "$D/j_panel_dmr$GEN/bl200_d10_cellline" ;;
    08)  echo "$D/j_panel_dmr${GEN}_panels" ;;
    09|10) [ -n "$ONE" ] && echo "$D/j_panel_dmr${GEN}_panels" || : ;;
    11)  echo "$D/j_pool${GEN}_cellline" ;;
    *)   : ;;
  esac
}

step(){ # step <번호> <제목> <명령...>
  local no=$1 title=$2; shift 2
  if [ -n "$ONE" ] && [ "$ONE" != "$no" ]; then return 0; fi
  say "── $no  $title"
  "$@" || { say "!! 실패 $no"; exit 1; }
}

# 검사는 「돌기 전에 한 번」 한다. 단계마다 검사하면 04 가 j_jsd 에 쓴
#   직후 04b 가 그 폴더를 「비어 있지 않다」고 막아 전체 실행이 항상 04b 에서 죽는다
#   (04·04b·04c 는 같은 폴더에 쓴다. 02 도 01 의 폴더에, 09·10 도 08 의 폴더에 쓴다).
#   그래서 「이번에 돌 단계」 가 쓰는 폴더를 모아 중복을 빼고 미리 한 번만 본다.
#   돌 단계 목록은 이 파일에서 직접 뽑는다. 손으로 적으면 단계를 더할 때 어긋난다.
_precheck() {
  local _ids _id _o _seen=" "
  if [ -n "$ONE" ]; then
    _ids=$ONE
  else
    # grep 이 실패하면 「검사한 폴더 0개」 로 조용히 통과했다
    #   (cat run_dmr.sh | bash 면 $0 가 bash 다). 비면 멈춘다 — 관문은 열린 채
    #   실패하면 안 된다.
    _ids=$(grep -oE "^step +[0-9a-z]+" "$0" 2>/dev/null | awk "{print \$2}")
    if [ -z "$_ids" ]; then
      echo "멈춤: 단계 목록을 못 읽었습니다 ($0). 파이프로 넘기지 말고" >&2
      echo "      bash $(basename "$0") 로 부르십시오 — 덮어쓰기 관문이 꺼집니다." >&2
      exit 1
    fi
  fi
  # 따옴표 없이 펼치면 METH_ROOT 에 공백이 있을 때 조각조각 나뉘어
   #   「없는 폴더」 로 전부 통과했다. 줄 단위로 읽는다.
  # 모르는 단계 번호(오타)면 검사도 실행도 없이 exit 0 이었다.
  #   「돌았는데 아무 일도 없음」은 성공으로 읽힌다. 먼저 막는다.
  if [ -n "$ONE" ]; then
    case " $(grep -oE "^step +[0-9a-z]+" "$0" 2>/dev/null | awk "{print \$2}" | tr "
" " ") " in
      *" $ONE "*) : ;;
      *) echo "멈춤: 그런 단계가 없습니다 — $ONE" >&2
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

say "DMR 시작: METH_ROOT=$METH_ROOT · GEN=$GEN"
_precheck

step 01 "후보 블록"        python -u "$S/01_candidates.py" --nd "$ND" --out "$D/j_candidates$GEN"
step 02 "블록별 CpG"       python -u "$S/02_blockcpg.py"  --nd "$ND" --cand "$D/j_candidates$GEN"
# JSD_LIMIT 을 주면 BAM 에서 읽을 리드 수를 제한한다 — 「도는지」 만 볼 때 쓴다.
#   본 실행에서는 주지 않는다. 전수를 읽어야 값이 맞는다.
_LIM=""
[ -n "${JSD_LIMIT:-}" ] && _LIM="--limit $JSD_LIMIT"
step 03 "창별 무늬수"      python -u "$S/03_jsdcount.py"  --nd "$ND" --cand "$D/j_candidates$GEN" \
                                --out "$D/j_jsdcount$GEN" $_LIM
step 04 "JSD"              python -u "$S/04_jsd.py"       --dir "$D/j_jsdcount$GEN" \
                                --cand "$D/j_candidates$GEN" --out "$D/j_jsd$GEN" --k "${PATTERN_K:-3}" --minr "${MIN_READS_PER_WINDOW:-6}"
# 04b·04c 를 사슬에 되돌렸다. 「감도 확인용」으로 보고 뺐는데 틀렸다 —
#   04b 가 null_*.parquet, 04c 가 boot_*.parquet 을 만들고, 06b 가 그 둘로
#   jsdb200 판(보고하는 세 판 중 하나)을 뽑는다. 빼면 08 이 jsdb200 에서 멈춘다.
step 04b "JSD 귀무 분포"    python -u "$S/04b_jsdnull.py"  --dir "$D/j_jsdcount$GEN" \
                                --cand "$D/j_candidates$GEN" --out "$D/j_jsd$GEN" --k "${PATTERN_K:-3}" --minr "${MIN_READS_PER_WINDOW:-6}"
step 04c "JSD 부트스트랩"   python -u "$S/04c_jsdboot.py"  --dir "$D/j_jsdcount$GEN" \
                                --cand "$D/j_candidates$GEN" --out "$D/j_jsd$GEN" --k "${PATTERN_K:-3}" --minr "${MIN_READS_PER_WINDOW:-6}"
step 05 "LLR 참조표"       python -u "$S/05_llrref.py"    --dir "$D/j_jsdcount$GEN" \
                                --cand "$D/j_candidates$GEN" --out "$D/j_llrref$GEN" --k "${PATTERN_K:-3}" --minr "${MIN_READS_PER_WINDOW:-6}"

_JF=$(ls "$D/j_jsd$GEN"/jsd_K${PATTERN_K:-3}_minr${MIN_READS_PER_WINDOW:-6}.parquet 2>/dev/null | head -1)
step 06 "JSD 상위 200 선정" python -u "$S/06_select.py" --cand "$D/j_candidates$GEN" \
                                --jsd "$_JF" --out "$D/j_panel$GEN" --n 200
step 06b "보완판 JSD 200 (jsdb200)" python -u "$S/06b_select_supplement.py"
step 07 "Baseline 순위 200" python -u "$S/07_panel.py"  --cand "$D/j_candidates$GEN" \
                                --out "$D/j_panel_dmr$GEN/bl200_d10_cellline" --topn 200

step 08 "검증용 형식으로 내보내기" python -u "$S/08_export_panels.py"
step 09 "귀무A 3판 (씨앗 20260914+i)"      python -u "$S/09_null_a.py"
step 10 "귀무B 3판 (씨앗 20260914+100+i)" python -u "$S/10_null_b.py"

# 11단계를 되돌렸다. 「하류가 안 읽는다」고 보고 뺐는데 틀렸다 —
#   11 이 만드는 j_pool<GEN>_cellline/pool_{gbm,normal}_all.bam 을
#   2_train/steps/01_...Mixing.py:56-57 이 그대로 읽는다(섞기의 원천 풀이다).
#   빼면 2단계 첫 걸음이 pysam FileNotFoundError 로 죽는다.
#   주의: run_pool.sh 는 panel_union3_<무엇>.csv 를 입력으로 받는데 그것을 만드는
#   코드는 여기 없다. 06_select 의 산출을 손으로 합친 것이다. README 에 적어 두었다.
step 11 "패널 영역 풀 BAM" bash "$S/run_pool.sh" cellline

say "DMR 끝: $D/j_panel_dmr${GEN}_panels/<패널>/DMR_confirmed_<패널>.csv"
