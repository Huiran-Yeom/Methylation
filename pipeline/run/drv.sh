#!/bin/bash
# 2·3·4단계 드라이버 하나(2026-10-27 정리판) — 옛 r20_run2.sh · r1_run2.sh · r2_run2.sh 와 보조 · 마무리 스크립트를 합쳤다.
#   쓰는 법: bash drv.sh <판> <코드 폴더> <main|help|finish> [기다릴 PGID …]
#     판     r20_gbm · r20_lung · r1_10 · r1_5 · r2_10 · r2_5 (profiles.sh)
#     main   대기열을 새로 만들고 → 자리 NL 개로 학습(T) · 검증(V) → 관문 → 4단계 표
#     help   돌고 있는 대기열에서 다음 항목을 꺼내 돈다(HID=<번호>). 대기열이 비면 끝난다. 관문 · 표는 안 한다.
#     finish 인자로 준 PGID 가 모두 끝나기를 기다린 뒤 관문 → 4단계 표(주 자리가 보조보다 먼저 끝나 «검증 미완» 으로 멈췄을 때)
#   환경 변수: DRY=1(무엇을 할지만 적음) · NV_SET · NL_SET · MEMW_SET · O_SET(대기열 폴더를 바꿈 · 시험용) · REPORT_ONLY_OUT(표를 이 폴더에 · 시험용) · FINISH_TEST=1(finish 를 관문까지만)
#   관문의 예외는 판 설정에 사전등록된 것(EXC_*)만 쓴다.
set -u
NAME=${1:?판 이름}; CODEDIR=${2:?코드 폴더 이름}; MODE=${3:?main · help · finish}; shift 3
for k in NL MEMW O; do         # 옛 드라이버의 이름(NL=1 bash …)은 판 설정이 덮어 조용히 무시된다 — 받지 않고 멈춘다
  [ -n "${!k:-}" ] && { echo "$k 는 받지 않는다 — ${k}_SET 를 쓴다"; exit 1; }
done
H=$(cd "$(dirname "$0")" && pwd)
. "$H/lib.sh"; . "$H/profiles.sh"
profile "$NAME" || exit 1
P=$HOME/pub/$CODEDIR/pipeline; export METH_CONF_DIR=$P/0_setup        # lib 뒤에 정한다(lib 은 P 를 건드리지 않는다)
[ -d "$P" ] || { echo "코드 폴더가 없다: $P"; exit 1; }
DRY=${DRY:-0}; HID=${HID:-}
R=$METH_ROOT
O=${O_SET:-$O}; mkdir -p "$O"; L=$O/run_$DS.log
NL=${NL_SET:-$NL}; NV=${NV_SET:-$NV}; MEMW=${MEMW_SET:-$MEMW}
Q=$O/q_$DS.txt; QL=$O/q_$DS.lock; RUN=$O/running_$DS; mkdir -p "$RUN"
TAG=""; [ "$MODE" = help ] && TAG="(보조 ${HID:-?}) "; [ "$MODE" = finish ] && TAG="(마무리) "
say(){ echo "[$(date '+%m-%d %H:%M:%S')] $TAG$*" >> "$L"; }
avail(){ df -BG --output=avail /ssd_data | tail -1 | tr -dc 0-9; }
mem_used(){ free -g | awk 'NR==2{print $3}'; }
mem_wait(){ while [ "$(mem_used)" -ge "$MEMW" ]; do sleep 60; done; }
del_ok(){ local r; for r in $DEL_ROOTS; do [ "$R" = "$r" ] && return 0; done; return 1; }
pop(){ ( flock 9; head -1 "$Q"; sed -i '1d' "$Q" ) 9> "$QL"; }
push_front(){ ( flock 9; { printf '%s\n' "$@"; cat "$Q"; } > "$Q.new" && mv "$Q.new" "$Q" ) 9> "$QL"; }

disk_guard(){
  [ "$(avail)" -ge 150 ] && return 0
  say "  디스크 $(avail)G < 150G — 채점 끝난 검증 칸의 섞은 BAM 만 지운다"
  local l s c v x
  if del_ok; then
    for l in "$TMP_ROOT"/val_*.log; do
      tail -n 1 "$l" | grep -q "끝 rc=0" || continue
      s=$(grep -o "채점된 검체 [0-9]*/[0-9]*" "$l" | tail -1 | grep -oE "[0-9]+/[0-9]+")
      [ -n "$s" ] && [ "${s%/*}" = "${s#*/}" ] || continue
      c=$(head -1 "$l" | awk '{print $3}'); v=$(head -1 "$l" | awk '{print $5}')
      for x in "$R"/${c}_val/*/results/$v/step03_Preprocessing_for_ML; do [ -d "$x" ] && [ ! -L "$x" ] && rm -rf -- "$x"; done
    done
  fi
  say "  지운 뒤 $(avail)G"
  [ "$(avail)" -ge 100 ]
}
do_train(){ local p=$1 d=$2 r=$3 v=j${GEN}$1_$2_uni l=$TMP_ROOT/train_$1_$2.log ns=4 c
  if [ "$DRY" = 1 ]; then say "  (DRY) 학습 $p $d (조각 $ns) · $P/2_train"; echo "$v" >> "$O/feasible_$DS.txt"
    local items=(); for c in $COHORTS; do items+=("V $c $v"); done; push_front "${items[@]}"; return 0; fi
  if ! { [ -s "$l" ] && tail -1 "$l" | grep -q "끝 rc=0"; }; then
    say "  학습 시작 $p $d (조각 $ns)"; NSLICE=$ns train_cell "$p" "$d" "$r"; say "  학습 끝 $p $d · $(tail -1 "$l" | cut -c1-40)"
  fi
  if [ "$(best_n "$R/results/panels/j${GEN}${p}/$v")" = "$EXPECT_BEST" ] && tail -1 "$l" | grep -q "끝 rc=0"; then
    echo "$v" >> "$O/feasible_$DS.txt"
    del_ok && { rm -rf -- "$R/sample_data/panels/j${GEN}${p}/$v"; say "  $v sample_data 중간 산출 지움"; }
    local items=(); for c in $COHORTS; do items+=("V $c $v"); done
    push_front "${items[@]}"      # 항목 하나 = 인자 하나 (따옴표 없이 넘기면 «V · 코호트 · 판» 이 세 줄로 쪼개진다)
  elif grep -q "풀 부족으로 중단" "$l" 2>/dev/null; then say "  $v 계산 불가(풀 부족 · 55-2)"; echo "$v" >> "$O/infeasible_$DS.txt"
  else say "  !! $v 학습 실패"; echo "$v" >> "$O/failed_$DS.txt"; fi
}
do_val(){ local c=$1 v=$2 l=$TMP_ROOT/val_$1_$2.log dd T0 rc
  [ "$DRY" = 1 ] && { say "  (DRY) 검증 $c $v · 검체 동시 $NV"; return 0; }
  [ -s "$l" ] && tail -1 "$l" | grep -q "끝 rc=0" && return 0
  disk_guard || { say "멈춤: 디스크"; return 1; }
  dd=${v%_uni}; dd=${dd##*_}
  ( export COHORT=$c VERSION=$v PANEL=$(echo "$v" | sed "s/^j${GEN}//;s/_.*//")
    case "$dd" in 5k) export TARGET_DEPTH_READS=5000;; 50k) export TARGET_DEPTH_READS=50000;; 100k) export TARGET_DEPTH_READS=100000;; esac
    echo "[$(date '+%m-%d %H:%M:%S')] $c · $v · METH_ROOT=$R" > "$l"
    T0=$(date +%s); bash "$P/3_validate/run_validate.sh" "$v" "$NV" >> "$l" 2>&1; rc=$?
    echo "[$(date '+%m-%d %H:%M:%S')] 끝 rc=$rc · $(( ($(date +%s)-T0)/60 ))분" >> "$l" )
  say "  검증 $c $v · $(grep -o '채점된 검체 [0-9]*/[0-9]*' "$l" | tail -1) · $(tail -1 "$l" | cut -c1-40)"
}
lane(){ local it
  while :; do
    it=$(pop)
    if [ -z "$it" ]; then
      [ -n "$(ls -A "$RUN")" ] || break          # 아무도 학습 중이 아니면 더 들어올 항목이 없다
      sleep 60; continue
    fi
    set -- $it
    [ "$#" -ge 3 ] || { say "  !! 대기열 항목이 이상하다 — 건너뜀: «$it»"; continue; }
    mem_wait
    case "$1" in
      T) touch "$RUN/$2_$3"; do_train "$2" "$3" "$4"; rm -f "$RUN/$2_$3" ;;
      V) do_val "$2" "$3" ;;
    esac
  done; }

exc_ok(){   # 사전등록 예외(EXC_*) 한 칸이 정말 그 칸인가 — 최신 채점 로그의 실패 검체 · 검체 로그 마지막 실행 구간의 «풀 부족»
  local pn=${EXC_V#j${GEN}}; pn=${pn%%_*}; local sl=$TMP_ROOT/logs/$EXC_C/j${GEN}$pn/$EXC_V
  local vl=$TMP_ROOT/val_${EXC_C}_${EXC_V}.log f n got t0
  t0=$(date -d "$(date +%Y)-$(head -1 "$vl" | cut -c2-15)" +%s 2>/dev/null) || { say "  사전등록 예외 아님: 검증 로그 시작 시각을 못 읽음"; return 1; }
  f=$(ls -t "$TMP_ROOT"/logs/$EXC_C/score_${EXC_V}_*.log 2>/dev/null | head -1); [ -n "$f" ] || return 1
  [ "$(stat -c %Y "$f")" -ge "$t0" ] || { say "  사전등록 예외 아님: 최신 채점 로그가 검증 시작 전 것"; return 1; }
  got=$(grep -o "^!! [A-Za-z0-9_]* 채점 실패" "$f" | awk '{print $2}' | sort -u | tr '\n' ' ' | sed 's/ $//')
  [ "$got" = "$EXC_N" ] || { say "  사전등록 예외 아님: 채점 실패 검체 «$got»"; return 1; }
  for n in $EXC_N; do
    [ -f "$sl/$n.log" ] && [ "$(stat -c %Y "$sl/$n.log")" -ge "$t0" ] || { say "  사전등록 예외 아님: $n 검체 로그가 없거나 검증 시작 전 것"; return 1; }
    awk '/^=== 시작 /{s=NR} {a[NR]=$0} END{if (!s) exit 1; for (i = s; i <= NR; i++) print a[i]}' "$sl/$n.log" | grep -q "풀 부족" \
      || { say "  사전등록 예외 아님: $n 마지막 실행 구간에 «풀 부족» 없음"; return 1; }
  done
}
gate(){   # 관문: 학습 실패 없음 · feasible 의 모든 칸이 n/n · 끝 rc=0 (사전등록 예외만 따로)
  [ -s "$O/failed_$DS.txt" ] && { say "멈춤: 학습 실패 칸 $(tr '\n' ' ' < "$O/failed_$DS.txt")"; return 1; }
  local bad=0 exc=0 c v l s
  for c in $COHORTS; do while read -r v; do l=$TMP_ROOT/val_${c}_${v}.log
    s=$(grep -o "채점된 검체 [0-9]*/[0-9]*" "$l" 2>/dev/null | tail -1 | grep -oE "[0-9]+/[0-9]+")
    tail -1 "$l" 2>/dev/null | grep -q "끝 rc=0" && [ -n "$s" ] && [ "${s%/*}" = "${s#*/}" ] && continue
    if [ -n "${EXC_C:-}" ] && [ "$c $v ${s:-}" = "$EXC_C $EXC_V $EXC_S" ] && exc_ok; then
      exc=$((exc + 1)); say "  사전등록 예외: $c $v $s · 풀 부족($EXC_N) — 1차 지표(5k) 밖 · 그대로 싣는다"; continue
    fi
    bad=1; say "  검증 미완 $c $v ${s:-}"
  done < <(sort -u "$O/feasible_$DS.txt"); done
  [ "$bad" = 0 ] || { say "멈춤: 검증 미완 칸"; return 1; }
  say "학습 · 검증 관문 통과$([ "$exc" -gt 0 ] && echo "(사전등록 예외 ${exc}칸)")"
}
report(){   # 4단계 — 원본(코호트 표 · cv_all.csv)을 먼저 끝낸다. 표는 임시 파일에 만든 뒤 rc 0 일 때만 바꿔 넣는다.
  #   Markdown 은 보기용이라 맨 마지막에 만든다. 실패해도 원본은 이미 있다.
  local out=${REPORT_ONLY_OUT:-$OUT} c pr rc n mdbad=0
  mkdir -p "$out/training" "$out/validation"
  for c in $COHORTS; do
    pr=$(person_re "$c")
    ( cd "$P/4_report" && if [ -n "$pr" ]; then COHORT=$c DEPTH=5k PERSON_RE="$pr" bash make_report.sh; else COHORT=$c DEPTH=5k bash make_report.sh; fi )       > "$out/validation/${c}_tables.txt.new" 2>&1; rc=$?
    [ "$rc" = 0 ] || { say "멈춤: 표 $c rc=$rc · $out/validation/${c}_tables.txt.new 를 본다"; return 1; }
    mv "$out/validation/${c}_tables.txt.new" "$out/validation/${c}_tables.txt"
    say "  표 $c rc=0 · $(wc -l < "$out/validation/${c}_tables.txt")줄"
  done
  python "$H/cv_all.py" "$R" "$GEN" "$out/training/cv_all.csv.new"; rc=$?
  [ "$rc" = 0 ] || { say "멈춤: 학습 표 rc=$rc"; return 1; }
  mv "$out/training/cv_all.csv.new" "$out/training/cv_all.csv"
  say "  학습 표 rc=0 · $out/training/cv_all.csv"
  n=$(grep -h "^수렴 경고 [0-9]*건" "$TMP_ROOT"/train_*.log 2>/dev/null | awk '{s+=$3+0}END{print s+0}'); say "  수렴 경고 합 ${n}건 (보고만 · CONVERGENCE_GATE 기본값)"
  say "  P_N 빈 창: $(grep -h "P_N : 정상 .*빈 창" "$TMP_ROOT"/train_*.log 2>/dev/null | sort -u | tr '
' ' ' | cut -c1-300)"
  for c in $COHORTS; do
    python "$H/tables_md.py" "$out/validation/${c}_tables.txt" > "$out/validation/${c}_tables.md" || { mdbad=1; say "  !! Markdown 실패 · 표 $c (원본 .txt 는 씀)"; }
  done
  python "$H/cv_md.py" "$out/training/cv_all.csv" > "$out/training/cv_all.md" || { mdbad=1; say "  !! Markdown 실패 · 학습 표 (원본 .csv 는 씀)"; }
  [ "$mdbad" = 0 ] || { say "멈춤: Markdown(원본은 씀)"; return 1; }
  say "  Markdown 씀(검증 표 · 학습 표)"
}

case "$MODE" in
  main|help)
    [ -s "$O/stage1_ok_$DS" ] || { say "멈춤: 관문 ① 표지 없음 — $STAGE1 부터"; exit 1; } ;;
esac
help_alive(){ local f p n=0; for f in "$O"/help_$DS.*.pid; do [ -e "$f" ] || continue; p=$(cat "$f"); [ -d "/proc/$p" ] && n=$((n + 1)) || rm -f "$f"; done; echo "$n"; }
case "$MODE" in
  main)   # 같은 판의 main 은 하나만 · 살아 있는 보조 자리가 있으면 대기열 · 표지를 지우지 않게 멈춘다
    exec 8> "$O/main_$DS.lock"; flock -n 8 || { say "멈춤: 같은 판의 main 이 이미 돈다"; exit 1; }
    # 보조 자리 등록과 직렬화한다 — 확인부터 대기열을 다 만들 때까지 등록 잠금을 쥔다
    exec 7> "$O/reg_$DS.lock"; flock 7
    [ "$(help_alive)" = 0 ] || { say "멈춤: 살아 있는 보조 자리 $(help_alive)개 — 끝난 뒤 main 을 띄운다"; exit 1; } ;;
  help)
    # main 이 대기열을 다시 만드는 중이면 등록 잠금에서 기다린다. 표지 이름은 PID 로 고유하게
    exec 7> "$O/reg_$DS.lock"; flock 7
    echo $$ > "$O/help_$DS.$$.pid"; trap 'rm -f "$O/help_$DS.$$.pid"' EXIT
    flock -u 7; exec 7>&- ;;
esac
case "$MODE" in
  main)
    say "── 2·3단계 · 코드 $P · 동시 $NL · 03 조각 $NSLICE · 검체 $NV · 디스크 $(avail)G · DRY=$DRY"
    ( exec 7>&- 8>&-; while kill -0 $$ 2>/dev/null; do u=$(mem_used); [ "$u" -ge 40 ] && echo "[$(date '+%m-%d %H:%M:%S')] !! 메모리 ${u}G" >> "$L"; sleep 120; done ) &
    # 대기열 · 상태 파일은 매번 새로 만든다 — 끝난 칸은 V 를 앞에, 안 끝난 칸은 T 를 뒤에, 풀 부족 칸은 infeasible 로.
    : > "$Q"; : > "$O/feasible_$DS.txt"; : > "$O/infeasible_$DS.txt"; : > "$O/failed_$DS.txt"
    rm -f "$RUN"/* 2>/dev/null
    [ "$DRY" = 1 ] || rm -f "$TMP_ROOT/val__.log" "$TMP_ROOT"/val_*__*.log 2>/dev/null   # 빈 이름 검증 로그(옛 판과 같음 · DRY 는 안 지움)
    CELLS=$( for d in "100k 100000" "50k 50000" "5k 5000"; do for p in $REAL; do echo "$p $d"; done; done
             for p in rand1200 rand2200 rand3200 randb1200 randb2200 randb3200; do echo "$p 5k 5000"; done )
    while read -r p d r; do
      v=j${GEN}${p}_${d}_uni; l=$TMP_ROOT/train_${p}_${d}.log
      if [ -s "$l" ] && tail -1 "$l" | grep -q "끝 rc=0" && [ "$(best_n "$R/results/panels/j${GEN}${p}/$v")" = "$EXPECT_BEST" ]; then
        echo "$v" >> "$O/feasible_$DS.txt"; for c in $COHORTS; do echo "V $c $v"; done >> "$Q"
      elif grep -q "풀 부족으로 중단" "$l" 2>/dev/null; then echo "$v" >> "$O/infeasible_$DS.txt"
      else echo "T $p $d $r" >> "$Q.t"; fi
    done <<< "$CELLS"
    [ -s "$Q.t" ] && cat "$Q.t" >> "$Q"; rm -f "$Q.t"
    flock -u 7; exec 7>&-                      # 대기열이 다 됐다 — 이제 보조 자리가 붙어도 된다
    say "대기열 $(wc -l < "$Q")항목 · 자리 $NL"
    pids=(); for i in $(seq 1 "$NL"); do lane "$i" & pids+=($!); sleep 5; done
    for i in "${pids[@]}"; do wait "$i"; done
    [ "$DRY" = 1 ] && { say "(DRY) 대기열 비움 · 판 $(sort -u "$O/feasible_$DS.txt" | wc -l) · 끝"; exit 0; }
    gate || exit 1
    report || exit 1 ;;
  help)
    say "── 보조 자리 · 코드 $P · 03 조각 $NSLICE · 검체 $NV · 디스크 $(avail)G · DRY=$DRY"
    say "대기열 $(wc -l < "$Q")항목 남음 · 보조 자리 시작"
    lane "h$HID"
    say "대기열 비음 · 보조 자리 끝" ;;
  finish)
    say "── 마무리 · 그룹 ${*:-없음} 이 끝나기를 기다린다"
    for g in "$@"; do while [ -n "$(ps -o pid= -g "$g" 2>/dev/null)" ]; do sleep 60; done; done
    say "── 마무리 · 그룹 모두 끝 — 관문 · 4단계"
    [ "$DRY" = 1 ] && { say "(DRY) 관문 · 표 · 끝"; exit 0; }
    gate || exit 1
    [ "${FINISH_TEST:-0}" = 1 ] && { say "(시험) 관문까지만"; exit 0; }
    report || exit 1 ;;
  *) echo "main · help · finish"; exit 1 ;;
esac
say "끝 · 디스크 $(avail)G"
