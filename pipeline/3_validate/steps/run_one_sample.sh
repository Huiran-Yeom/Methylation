#!/bin/bash
# 설정은 0_setup/config.conf 에서만 고친다.
#   2026-09-25 식별자를 ASCII 로 바꿨다. 이 셸은 한글 변수명에서 죽는다
#   (`_여기=...: No such file or directory`). 파일 이름도 경로.conf -> config.conf.
_HERE=$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)
# 환경변수로 주는 값이 이긴다. meth_config.py 의 값() 과 같은 순서로 맞춘다.
#   conf 를 그냥 source 하면 환경변수를 덮어써서, 환경변수로 격리한 실행이
#   조용히 실 경로로 넘어간다(09-26 에 실제로 났다). 미리 보관하고 되돌린다.
#   값은 eval 에 태우지 않는다. 이름만 고정 목록으로 돈다.
# 2026-09-28: 보관 대상을 손으로 들고 있으면 목록에 없는 키는 conf 가 환경변수를
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
# 2026-09-25: 못 찾았을 때 조용히 넘어가면 VERSION·PANEL·SEED 가 아예 안 잡히고,
#   환경에 남아 있던 METH_CONF_DIR 로 「다른 트리의 설정」 을 읽는다. 죽게 한다.
: "${METH_CONF_DIR:?config.conf 를 못 찾았다 — 트리 밖에서 돌리려면 METH_CONF_DIR 를 직접 주라}"
: "${GEN:?GEN(패널 세대)이 없다 — config.conf 의 GEN 을 채우십시오}"
: "${METH_ROOT:=/ssd_data/Methylation}"
: "${DATASET_ROOT:=/ssd_data/dataset}"
# 2026-09-29: TMP_ROOT 를 여기서도 지킨다. 드라이버(run_validate.sh)만
#   지키고 로그를 「쓰는」 이쪽이 ~ 를 박아 두어, 그 값을 바꾸면 실패 안내가
#   영원히 비어 있을 폴더를 가리켰다.
: "${TMP_ROOT:=$HOME/tmp}"
# 2026-09-29: 기본값을 저장소 안으로 돌린다. 예전 기본값(~/Methylation_rev1/...)은
#   받아서 돌리는 사람에게 없는 경로라 학습 스크립트를 못 찾은 채 돌았다.
: "${TRAIN_SRC:=$_HERE/../../2_train/steps}"
: "${ANALYSIS_DIR:=$_HERE}"

# j 판 외부검증: 검체 하나를 끝까지 돌린다. 2026-08-24
#   사용:  run_one_sample.sh <검체이름> <판이름>
#   폴더가 없으면 make_val_config.py 로 만든다. 각 단계는 산출물이 있으면 건너뛴다.
set -u
# ─── 2026-09-15 [비우기] 죽이지 말고 비운다 ───────────────────────
#   자원이 다툴 때 드라이버를 죽이면 CSV 를 쓰던 자식이 **잘린 파일**을 남기고,
#   그 잘린 파일이 존재 검사를 통과한다(오늘 실제로 만들 뻔했다).
#   대신 표시 파일 하나만 보고 **새로 걸지 않는다**. 도는 것은 끝까지 간다.
#
#     멈추기  touch ~/tmp/lock/STOP
#     풀기    rm   ~/tmp/lock/STOP
#
#   검체마다 새 셸이 뜨므로 여기 한 곳이면 모든 런처(이어달리기·보조 xargs)를 덮는다.
# 2026-09-17 [OOM 우선순위] 스왑이 0 이라(실측) 메모리가 마르면 곧장 OOM 이거나 fork 실패다.
#   권한 없는 프로세스는 자기 oom_score_adj 를 **올릴 수만 있고 내릴 수 없다** —
#   그래서 감시를 보호하는 대신 **워커를 먼저 죽게** 만든다. 자식(python·bismark)이 물려받는다.
#   09-17 14:0x 에 sshd 가 악수를 못 해 1시간 40분 접속 불가였다. 그때 죽었어야 할 것은 워커다.
echo 800 > /proc/self/oom_score_adj 2>/dev/null
if [ -f "$TMP_ROOT/lock/STOP" ]; then
  echo "[$(date +%H:%M:%S)] 그만 표시 — $1 $2 새로 걸지 않습니다" >> "$TMP_ROOT/drop_caches.log"
  exit 0
fi
# ──────────────────────────────────────────────────────────────
s="$1"; v="$2"
# 2026-09-29: 기본값을 없앴다. 우리 내부 코호트 이름이 기본값이면, 이 스크립트를
#   직접 부르는 사람(문서에 있는 사용법이다)이 COHORT 를 빠뜨렸을 때 조용히
#   「남의 코호트」 로 돌고 그 트리에 쓴다. 안 주면 멈춘다.
: "${COHORT:?COHORT 를 주십시오 — 3_validate/run_validate.sh 와 같은 값이어야 합니다}"
C=$COHORT
# 2026-09-29: WORK_ROOT 를 따른다. 예전에는 여기만 ~ 를 박아 두고 채점 쪽은
#   WORK_ROOT 를 봐서, 그 값을 바꾸면 피처와 채점이 서로 다른 자리를 봤다.
W=${WORK_ROOT:-$HOME}/${C}_val/$v/$s
# 2026-09-25: 코호트를 자리에 넣는다. 없으면 재현시험과 실운영이 같은 판·같은 검체를
#   돌 때 「한 로그」 에 겹쳐 쓴다. 실제로 겹쳐서 옛 실패 줄을 새 것으로 오독했다.
D=$TMP_ROOT/logs/${C}/${v%%_*}/${v}
# 옛 로그 자리(코호트 이름이 붙기 전)를 이어 쓰던 분기다. 코호트 이름을
#   환경변수로 받게 된 뒤로는 새 자리만 쓴다.
O=
# 08-26: 로그를 판별 폴더로 내린다.
#   옛 자리에 이미 있으면 그대로 쓴다. 세는 곳이 둘로 갈리지 않게.
if [ -n "$O" ] && [ -f "$O" ]; then L="$O"; else mkdir -p "$D"; L="$D/${s}.log"; fi
# ---- 2026-09-29 [잠금을 먼저] 검체+판 단위 잠금 ----
#   예전에는 config.py 생성이 잠금 「앞」 에 있었다. 드라이버 둘이 같은 검체를 잡으면
#   둘 다 make_val_config 로 들어가 같은 config.py 를 동시에 쓴다. 잘린 config.py 는
#   경로를 검증 트리로 돌리는 줄(OVERRIDE 끝)에 못 닿아 「학습 트리」 를 가리킨다 —
#   즉 검증이 학습 산출물을 덮는다. 잠금을 먼저 잡으면 그 경주가 없다.
#   실패하면 조용히 나간다. 다른 쪽이 하고 있으니 그게 맞다.
mkdir -p "$TMP_ROOT/lock/val" "$W"
exec 9>"$TMP_ROOT/lock/val/${v}__${s}.lock"
if ! flock -n 9; then
  echo "=== 다른 실행이 잡고 있어 건너뜁니다 $(date +%H:%M:%S) $s / $v" >> "$L"
  exit 0
fi
if [ ! -f "$W/config.py" ]; then
  # 2026-09-24 make_val_config.py 가 있는 곳. config.conf 의 ANALYSIS_DIR.
  cd "$ANALYSIS_DIR" || exit 1
  python make_val_config.py "$C" "$s" "$v" >> "$L" 2>&1 || { echo "!! config 실패 $s" >> "$L"; exit 1; }
fi
cd "$W" || exit 1

# ---- 2026-09-15 [이미 끝난 검체 건너뛰기] ----
#   한 검체가 5~10분이라 178검체를 다시 도는 것은 몇 시간이다.
#   **산출물로 판정한다**: 로그는 이어붙기라 옛 실행 것이 남아 있을 수 있다.
#   피처 5종이 다 있고 개수가 같으면 끝난 것이다. 다시 돌리려면 VALREDO=1.
if [ "${VALREDO:-0}" != "1" ]; then
  _RDONE=$METH_ROOT/${C}_val/$s/results/$v/step04_ML_classifier
  _dn=1; _dc=""
  for _k in mean entropy jsd pdr llr; do
    _m=$(ls $_RDONE/[0-9]*_${_k}_*.csv 2>/dev/null | wc -l)
    [ "$_m" -ge 1 ] || { _dn=0; break; }
    if [ -z "$_dc" ]; then _dc=$_m; elif [ "$_m" != "$_dc" ]; then _dn=0; break; fi
  done
  if [ "$_dn" = "1" ]; then
    echo "=== 건너뜀 $(date +%H:%M:%S) $s / $v — 피처 5종 x ${_dc}비율 이미 있음" >> "$L"
    exit 0
  fi
fi

# 2026-08-24: val04_1c 는 08-21 에 새로 만든 것이라 그 전에 생긴 폴더엔 없다.
#   없으면 복사한다. 원본은 판 폴더가 아니라 한 단계 위에 있다.
if [ ! -f val04_1c.py ]; then
  cp -f "$TRAIN_SRC/06_jsdfeat.py" val04_1c.py 2>/dev/null \
    || echo "!! val04_1c 원본 없음" >> "$L"
fi
# ---- 2026-09-15 [스크립트 최신화] 검체 폴더가 낡는 것을 막는다 ----
#   위 17줄은 config.py 가 있으면 make_val_config 를 다시 안 부른다. 그래서 먼저
#   만들어진 검체 폴더에 **옛 스크립트가 남는다.** 실제로 val03_2.py 가 옛 판본이라
#   _bgreads 가 없어 KeyError 로 죽었다(학습 쪽에서 이미 고친 결함이었다).
#   코드를 고칠 때마다 검체 폴더가 낡으므로 지나가는 자리에서 한 번 막는다.
#   config.py 는 검체마다 다르므로 건드리지 않는다. 스크립트만 맞춘다.
# 2026-09-24 학습 스크립트 원본. config.conf 의 TRAIN_SRC 로 바꾼다.
_SRCDIR=$(ls -d "$TRAIN_SRC"/*/$v 2>/dev/null | head -1)
[ -n "$_SRCDIR" ] || _SRCDIR=$(ls -d "$TRAIN_SRC" 2>/dev/null | head -1)
if [ -n "$_SRCDIR" ]; then
  # 2026-09-24 이름 짝을 config.conf 의 STEP_MAP 에서 읽는다.
  #   비어 있으면 아래 기본값을 쓴다. 파일명이 바뀌면 설정만 고친다.
  _pairs="${STEP_MAP:-val03_1:01_Preprocessing_for_ML_RawData_Mixing.py
val03_2:02_Preprocessing_for_ML_RawData_Mixing.py
val03_3:03_Preprocessing_for_ML_RawData_Mixing.py
val04_0:04_coverage_check.py
val04_1:05_Feature_Matrix.py
val04_1c:06_jsdfeat.py
val04_1d:07_llrfeat.py}"
  # 2026-09-25 헤어독을 버리고 파이프로 읽는다. 헤어독 안의 $( ) 치환이
  #   세미콜론 목록과 얽혀 깨졌다. 파이프는 한 줄이고 깨지지 않는다.
  _MISS=$(mktemp)
  printf '%s\n' "$_pairs" | tr ';' '\n' | while IFS=: read -r _n _f; do
    [ -n "$_n" ] || continue
    _src="$_SRCDIR/$_f"; [ -f "$_src" ] || _src="$(dirname "$_SRCDIR")/$_f"
    # 2026-09-28 [★ 조용한 건너뛰기] 이름이 안 맞으면 옛 val04_1x.py 가 그대로 돈다.
    #   그 피처만 빠진 채 완주하고 채점이 결측을 평균으로 메워 AUC 0.500 이 난다.
    #   while 이 파이프 서브셸이라 여기서 exit 해도 부모가 안 끝난다 — 표시를 남긴다.
    if [ ! -f "$_src" ]; then
      echo "$_f" >> "$_MISS"
      continue
    fi
    _a=$(md5sum "$_src" | cut -d" " -f1)
    _b=$(md5sum "$_n.py" 2>/dev/null | cut -d" " -f1)
    # 2026-09-29: cp 실패(디스크 참·권한)를 삼키면 「옛 스크립트가 그대로 돈다」.
    #   md5sum 이 없으면 _a·_b 가 둘 다 비어 「같다」 가 되어 아무것도 최신화 안 된다.
    if [ -z "$_a" ]; then echo "$_f (md5sum 실패)" >> "$_MISS"; continue; fi
    if [ "$_a" != "$_b" ]; then
      if cp -f "$_src" "$_n.py"; then
        echo "--- 최신화 $_n.py (${_b:0:8} -> ${_a:0:8})" >> "$L"
      else
        echo "$_f (복사 실패)" >> "$_MISS"
      fi
    fi
  done
  if [ -s "$_MISS" ]; then
    echo "!! 최신화할 원본이 없습니다 — STEP_MAP 의 이름과 원본 파일명이 어긋났습니다:" >> "$L"
    sed 's/^/     /' "$_MISS" >> "$L"
    echo "   그대로 두면 옛 스크립트가 돌아 그 피처만 빠진 채 완주합니다 (AUC 0.500)." >> "$L"
    rm -f "$_MISS"; exit 1
  fi
  rm -f "$_MISS"
else
  # 2026-09-29: 여기서 경고만 하고 계속하면 안 된다. 위 _MISS 분기와 「같은 실패」 다 —
  #   검체 폴더에 남아 있는 옛 스크립트가 그대로 돌아 그 피처만 빠진 채 완주하고,
  #   채점이 결측을 학습 평균으로 메워 AUC 0.500 이 난다. 그쪽은 exit 1 인데
  #   이쪽만 경고였다. 같은 결말에는 같은 처리를 한다.
  echo "!! 원본 폴더를 못 찾았습니다. TRAIN_SRC=$TRAIN_SRC" >> "$L"
  echo "   학습 스크립트를 최신화할 수 없습니다. 그대로 두면 옛 스크립트가 돌아" >> "$L"
  echo "   그 피처만 빠진 채 완주합니다 (AUC 0.500). config.conf 의 TRAIN_SRC 를 보십시오." >> "$L"
  exit 1
fi
echo "=== 시작 $(date +%H:%M:%S) $s / $v" >> "$L"
# 2026-08-24: val04_1c(JSD 피처)를 붙였다. 칸 목록은 학습이 고른 것을 그대로 쓴다.
#   검증 자료를 보고 다시 고르면 BEST.joblib 의 열과 안 맞는다.
export JSD_COLS=$METH_ROOT/results/panels/${v%%_*}/$v/step04_ML_classifier/jsd_cols_$v.txt
# 2026-09-29: 없으면 멈춘다. 예전에는 JSD_COLS="" 로 내려가 검증 검체의 결측을 보고
#   칸을 「다시」 골랐고, BEST.joblib 의 열과 어긋난 자리를 채점이 학습 평균으로 메워
#   AUC 0.500 이 났다. 위 두 분기와 같은 결말이다. 이제 학습 06 이 이 파일을 만든다.
# 2026-09-29: -f 가 아니라 -s 다. 학습이 칸을 하나도 못 고르면(커버리지가 얕거나
#   MAXNA 를 잘못 주면) 0바이트 파일이 생기고, 그것이 관문을 통과해 피처가 빈 채
#   완주한다. 막으려던 AUC 0.500 그대로다.
if [ ! -s "$JSD_COLS" ]; then
  echo "!! 칸 목록이 없습니다. $JSD_COLS" >> "$L"
  echo "   학습(2_train/run_train.sh 의 06)이 이 파일을 만듭니다. 학습을 먼저 돌리십시오." >> "$L"
  echo "   없이 돌리면 검증 자료로 칸을 다시 골라 AUC 0.500 이 납니다." >> "$L"
  exit 1
fi
# 2026-09-14: val04_1d(LLR) 를 붙인다. 목록에 없어 **검증이 llr 을 아예 안 만들었다.**
if [ ! -f val04_1d.py ]; then
  cp -f "$TRAIN_SRC/07_llrfeat.py" val04_1d.py 2>/dev/null     || echo "!! val04_1d 원본 없음" >> "$L"
fi
# 참조표와 후보 풀은 **판 세대**를 따라간다. 기본값에 기대면 학습/적용이 어긋난다.
#   학습(돌리기_j15.sh)은 llrref15 + j_candidates15 를 쓴다. 검증도 같아야 한다.
# 2026-09-26: 'j15*' 를 못 박으면 GEN 을 바꿨을 때 이 분기가 「안 걸리고」
#   아래 기본값(세대 없는 폴더)으로 내려간다. GEN 을 따른다.
# 2026-09-29: case 를 없앤다. 예전에는 판 이름이 j<GEN> 으로 시작할 때만 세대를
#   따라가고, 그 밖의 이름은 「세대 없는」 폴더로 내려갔다. 그 대체 분기에서는
#   JSD_CAND·JSD_JCOUNT 를 아예 안 내보내서, 06_jsdfeat 이 자기 기본값(세대 없는
#   j_candidates)으로 조용히 돌 수 있었다. 09-23 사고와 같은 모양이다.
#   판 이름과 무관하게 GEN 을 따른다. 넷을 한 자리에서 내보내 빠뜨릴 수 없게 한다.
export LLR_REF=${LLR_REF:-$METH_ROOT/results/dmr/j_llrref${GEN}/llrref_K${PATTERN_K:-3}_minr${MIN_READS_PER_WINDOW:-6}.parquet}
export LLR_CAND=${LLR_CAND:-$METH_ROOT/results/dmr/j_candidates${GEN}}
# val04_1c(06_jsdfeat)는 인자를 안 받으므로 환경변수로만 줄 수 있다.
export JSD_CAND=${JSD_CAND:-$METH_ROOT/results/dmr/j_candidates${GEN}}
export JSD_JCOUNT=${JSD_JCOUNT:-$METH_ROOT/results/dmr/j_jsdcount${GEN}}
[ -f "$LLR_REF" ] || { echo "!! 참조표 없음 $LLR_REF" >> "$L"; exit 1; }
# 2026-09-29: 파일 이름을 박지 않는다 — 설정의 K·MINR 을 바꾸면 학습(run_train.sh)과
#   여기가 「다른 문턱의 참조표」 를 보게 되고, 두 파일은 상위폴더만 다르다.
# 세대마다 이름이 같으므로 basename 으로는 구분이 안 된다. 상위폴더까지 찍는다.
echo "--- 참조 $(basename $(dirname $LLR_REF))/$(basename $LLR_REF) · 후보 $(basename $LLR_CAND)" >> "$L"
# 2026-09-29: 도는 순서를 「박지 않는다」. 여기가 STEP_MAP 의 네 번째 사본이었다 —
#   설정에 단계를 더해도 복사만 되고 「돌지 않아」, 그 피처가 조용히 빠졌다.
#   위에서 이미 만든 $_pairs(= STEP_MAP 또는 그 기본값)의 왼쪽 이름을 순서대로 쓴다.
_STEPS=$(printf '%s\n' "$_pairs" | tr ';' '\n' | cut -d: -f1 | grep -v '^$' | tr '\n' ' ')
[ -n "$_STEPS" ] || { echo "!! 돌 단계 목록이 비었다 — STEP_MAP 을 보십시오" >> "$L"; exit 1; }
echo "--- 단계 $_STEPS" >> "$L"
for st in $_STEPS; do
  # 2026-09-29: 없으면 조용히 넘기지 않는다. 위 최신화 루프가 원본을 확인했으므로
  #   여기서 없는 것은 「복사가 안 된」 것이고, 넘기면 그 피처만 빠진 채 완주한다.
  [ -f "$st.py" ] || { echo "!! $st.py 가 검체 폴더에 없다 — 복사가 안 됐다" >> "$L"; exit 1; }
  # 2026-09-29: 메모리 문을 왼쪽 이름(val03_2)이 아니라 「원본 파일명」 에 건다.
   #   왼쪽 이름은 STEP_MAP 에서 바꿀 수 있는 값이고(config.conf 가 그렇게 안내한다),
   #   그 이름을 바꾸면 문이 조용히 사라진다. 메모리를 쓰는 것은 02_...Mixing.py 다.
  _SRCF=$(printf %s "$_pairs" | tr ';' '\n' | awk -F: -v k="$st" '$1==k{print $2}' | head -1)
  case "$_SRCF" in 02_*) _GATE=1 ;; *) _GATE=0 ;; esac
  # 2026-09-17 [메모리 문] val03_2(혼합)만 검체당 약 1.35G 를 쓴다(실측 14:36).
  #   나머지 단계는 가볍다. **이 단계만 막으면 전체 병렬을 올려도 안전하다.**
  #   09-17 에 66워커가 동시에 이 단계로 몰려 62G 를 다 먹고 sshd 까지 막혔다.
  #   부하·프로세스 수로는 안 걸렸다 - 그래서 감시가 아니라 **여기**에 문을 단다.
  if [ "$_GATE" = 1 ]; then
    # 2026-09-17 [문 토큰화] 옛 코드는 도는 val03_2 를 **세고** 25 미만이면 통과시켰다.
    #   세기는 통과와 등록 사이에 틈이 있어 대기자들이 동시에 깨서 동시에 들어간다
    #   (sleep 20 이 정렬돼 있어 더 나쁘다). 25 는 상한이 아니라 관측값이었다.
    #   → **세지 않고 잡는다.** 잠금 25개를 토큰으로 두고 하나를 쥐어야 들어간다.
    #     잡은 수가 곧 도는 수다. 죽으면 커널이 fd 를 닫아 자동 반납된다.
    #   ★ 스왑이 0 이고(실측) 부하감시는 이 스크립트를 죽이지 않는다(renice 만).
    #     그래서 이 문이 **유일한 브레이크**다. 25 를 올리지 말 것:
    #     25 x 1.35G = 33.8G + 다른 단계·OS 약 8G = 42G → 남는 20G 가 sshd·감시의 생존선이다.
    _got=0
    for _i in $(seq 1 540); do
      _a=$(free -g | awk '/^Mem:/{print $7}')
      case "$_a" in ''|*[!0-9]*) _a=0 ;; esac      # 값이 이상하면 「부족」 으로 본다
      if [ "$_a" -ge "${MEMGATE_FREE:-22}" ]; then
        for _k in $(seq 1 ${MEMGATE_N:-22}); do
          # ★ fd 9 는 이 스크립트가 **검체별 잠금**에 이미 쓴다(위 exec 9>).
          #   여기서 9 를 재할당하면 그 잠금이 풀려 두 드라이버가 같은 검체를 잡는다.
          #   그래서 문은 fd 6 을 쓴다.
          exec 6>>"$TMP_ROOT/lock/val03_2.$_k"
          flock -n 6 && { _got=1; break; }
        done
      fi
      [ "$_got" = 1 ] && break
      [ "$_i" = 1 ] && echo "--- val03_2 문 대기 (여유 ${_a}G)" >> "$L"
      sleep $((15 + RANDOM % 11))                  # 정렬을 깬다
    done
    if [ "$_got" != 1 ]; then
      # 그냥 통과시키면 대기자들이 한꺼번에 몰려 09-17 14:0x 고장이 재현된다.
      # 이 검체를 실패로 두고 나간다. 피처가 없으므로 「건너뛰기」 에 안 걸려 다음에 다시 잡힌다.
      echo "!! val03_2 문 3시간 초과: 이 검체를 건너뛴다 (나중에 다시 잡힌다)" >> "$L"
      exit 1
    fi
  fi
  echo "--- $st $(date +%H:%M:%S)" >> "$L"
  # 2026-09-29: fd 6(메모리 토큰)·9(검체 잠금)를 자식에게 물려주지 않는다.
  #   물려주면 고아 자식이 잠금을 쥔 채 남아 토큰이 영영 반납되지 않고,
  #   그 검체는 영원히 「다른 실행이 잡고 있다」 로 보인다.
  python -u "$st.py" 6>&- 9>&- >> "$L" 2>&1 || { echo "!! $st 실패" >> "$L"; exit 1; }
  # 문 토큰 반납. 안 하면 뒤 단계까지 쥐고 있어 처리량이 25 로 묶인다.
  [ "$_GATE" = 1 ] && exec 6>&-         # 문 토큰만 반납. 검체 잠금(fd 9)은 그대로 둔다
done
echo "=== 끝 $(date +%H:%M:%S)" >> "$L"

# ---- 2026-09-15 [중간물 정리] 확인한 뒤에만 지운다 ----
#   옛 코호트는 검체당 7.3G 였고 178검체면 1.3TB 다. 정리가 없으면 45검체에서 디스크가 찬다.
#   피처 CSV 7종이 다 있고 개수가 서로 같을 때만 지운다. 없는데 지우면 Bismark 부터 다시다.
#   results/ 는 안 건드린다. sample_data/panels/<앞머리>/<판>/ 만 지운다.
#   끄려면 VALKEEP=1 로 부른다.
if [ "${VALKEEP:-0}" != "1" ]; then
  _SD=$METH_ROOT/${C}_val/$s/sample_data/panels/${v%%_*}/$v
  _RD=$METH_ROOT/${C}_val/$s/results/$v/step04_ML_classifier
  _ok=1; _cnt=""
  for _k in mean entropy jsd pdr llr; do
    _n=$(ls $_RD/[0-9]*_${_k}_*.csv 2>/dev/null | wc -l)
    [ "$_n" -ge 1 ] || { echo "   정리 안 함 — $_k 피처가 없습니다" >> "$L"; _ok=0; break; }
    if [ -z "$_cnt" ]; then _cnt=$_n
    elif [ "$_n" != "$_cnt" ]; then
      echo "   정리 안 함. 피처 개수가 다릅니다 ($_k $_n vs $_cnt)" >> "$L"; _ok=0; break
    fi
  done
  if [ "$_ok" = "1" ] && [ -d "$_SD" ]; then
    _sz=$(du -sm "$_SD" 2>/dev/null | cut -f1)
    rm -rf "$_SD"
    echo "   중간물 정리 $_SD (${_sz}M · 피처 5종 x ${_cnt}비율 확인)" >> "$L"
  fi
  # ---- 2026-09-17 [BAM 회수] 디스크가 쪼들릴 때만 지운다 (사용자 지시: "부족하면 지워") ----
  #   results/<판>/step03_.../mut_*/sampled_reads_*.bam 이 영구 산출물의 거의 전부다.
  #   50k 약 398M/칸 · 100k 약 696M/칸. 남은 1,668칸이면 912G 인데 여유가 1,008G 였다.
  #   [확인] 이 BAM 을 읽는 것은 val03_2·val03_3·val04_1·val04_1b/1c/1d 뿐이다.
  #          채점(score_one_sample.py)은 CSV 만 읽는다. 실측으로 확인했다.
  #          그러므로 피처 5종이 다 나온 뒤에는 아무도 안 읽는다.
  #   [잃는 것] 피처 계산만 다시 돌릴 수 없다. bismark 부터 다시다(검체당 약 40분).
  #          그래서 여유가 있을 땐 남긴다. DISKMIN 아래일 때만 회수한다.
  #   [남기는 것] background_map_*.csv · txt.tar.gz — "몇 명을 섞었나" 의 증거다. 09-17 에 "3명" 을 이걸로 잡았다.
  if [ "$_ok" = "1" ]; then
    _free=$(df -BG --output=avail "$METH_ROOT" 2>/dev/null | tail -1 | tr -dc '0-9')
    case "$_free" in ''|*[!0-9]*) _free=99999 ;; esac
    if [ "$_free" -lt "${DISKMIN:-200}" ]; then
      _BD=$METH_ROOT/${C}_val/$s/results/$v/step03_Preprocessing_for_ML
      _bn=$(ls "$_BD"/mut_*_reads/sampled_reads_*.bam 2>/dev/null | wc -l)
      if [ "$_bn" -gt 0 ]; then
        _bsz=$(du -cm "$_BD"/mut_*_reads/sampled_reads_*.bam 2>/dev/null | tail -1 | cut -f1)
        rm -f "$_BD"/mut_*_reads/sampled_reads_*.bam
        echo "   BAM 회수 ${_bn}개 ${_bsz}M (여유 ${_free}G < ${DISKMIN:-200}G)" >> "$L"
      fi
    fi
  fi
fi
