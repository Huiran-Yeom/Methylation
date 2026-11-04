#!/bin/bash
# 설정은 0_setup/config.conf 에서만 수정
# 변수명, 파일 이름은 ASCII 만 (한글 변수명이면 셸이 죽음)
_HERE=$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)
# config.conf 읽기. 환경변수가 우선 (meth_config.py 와 같은 규칙)
# conf 키를 전부 미리 보관해 뒀다가 conf 읽은 뒤 되돌림
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
# config.conf 못 찾으면 중단
: "${METH_CONF_DIR:?config.conf 를 못 찾았다 — 트리 밖에서 돌리려면 METH_CONF_DIR 를 직접 주라}"
# BLAS 스레드 고정 - 로지스틱 회귀 값이 스레드 수 따라 조금 달라짐
# 환경변수에 이미 있으면 그 값 사용, BLAS_THREADS=none 이면 고정 안 함
_bt=${BLAS_THREADS:-1}
case "$_bt" in
  [Nn][Oo][Nn][Ee]) ;;
  ''|*[!0-9]*|0) echo "멈춤: BLAS_THREADS 는 none 또는 1 이상의 정수여야 합니다: $_bt" >&2; exit 1 ;;
  *) export OMP_NUM_THREADS=${OMP_NUM_THREADS:-$_bt} OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-$_bt} MKL_NUM_THREADS=${MKL_NUM_THREADS:-$_bt}
     [ "$OMP_NUM_THREADS" = "$OPENBLAS_NUM_THREADS" ] && [ "$OMP_NUM_THREADS" = "$MKL_NUM_THREADS" ]        || echo "주의: 스레드 값이 갈린다 · OMP=$OMP_NUM_THREADS OPENBLAS=$OPENBLAS_NUM_THREADS MKL=$MKL_NUM_THREADS" >&2 ;;
esac
: "${GEN:?GEN(패널 세대)이 없다 — config.conf 의 GEN 을 채우십시오}"
: "${METH_ROOT:=/ssd_data/Methylation}"
: "${DATASET_ROOT:=/ssd_data/dataset}"
# TMP_ROOT 는 드라이버(run_validate.sh)와 같은 값이어야 함
: "${TMP_ROOT:=$HOME/tmp}"
# 기본값은 저장소 안 학습 스크립트
: "${TRAIN_SRC:=$_HERE/../../2_train/steps}"
: "${ANALYSIS_DIR:=$_HERE}"

# 외부검증: 검체 하나를 끝까지 돌림
#   사용법:  run_one_sample.sh <검체이름> <판이름>
#   폴더 없으면 make_val_config.py 로 만듦. 끝난 검체는 건너뜀
set -u
# 멈추기: 죽이지 말고 STOP 표시 - 새 검체만 안 걸고 도는 건 끝까지 감
#   (죽이면 쓰던 CSV 가 잘린 채 남음)
#     멈추기  touch ~/tmp/lock/STOP
#     풀기    rm   ~/tmp/lock/STOP
# OOM 때 워커가 먼저 죽게 oom_score_adj 올림 (스왑 0, 자식도 물려받음)
echo 800 > /proc/self/oom_score_adj 2>/dev/null
if [ -f "$TMP_ROOT/lock/STOP" ]; then
  echo "[$(date +%H:%M:%S)] 그만 표시 — $1 $2 새로 걸지 않습니다" >> "$TMP_ROOT/drop_caches.log"
  exit 0
fi
s="$1"; v="$2"
# COHORT 기본값 없음, 없으면 중단
: "${COHORT:?COHORT 를 주십시오 — 3_validate/run_validate.sh 와 같은 값이어야 합니다}"
C=$COHORT
# 검체 작업 폴더 (채점 루프와 같은 WORK_ROOT)
W=${WORK_ROOT:-$HOME}/${C}_val/$v/$s
# 로그 폴더. 코호트, 판별로 나눔
D=$TMP_ROOT/logs/${C}/${v%%_*}/${v}
# 옛 로그 자리 (비워 둠, 새 자리만 사용)
O=
if [ -n "$O" ] && [ -f "$O" ]; then L="$O"; else mkdir -p "$D"; L="$D/${s}.log"; fi
# 검체+판 단위 잠금. config.py 만들기 전에 잡음
#   (둘이 동시에 쓰다 잘린 config.py 는 학습 트리를 가리켜 학습 산출물을 덮음)
#   못 잡으면 다른 실행이 도는 중이라 그냥 나감
mkdir -p "$TMP_ROOT/lock/val" "$W"
exec 9>"$TMP_ROOT/lock/val/${v}__${s}.lock"
if ! flock -n 9; then
  echo "=== 다른 실행이 잡고 있어 건너뜁니다 $(date +%H:%M:%S) $s / $v" >> "$L"
  exit 0
fi
if [ ! -f "$W/config.py" ]; then
  # make_val_config.py 자리 (config.conf 의 ANALYSIS_DIR)
  cd "$ANALYSIS_DIR" || exit 1
  python make_val_config.py "$C" "$s" "$v" >> "$L" 2>&1 || { echo "!! config 실패 $s" >> "$L"; exit 1; }
fi
cd "$W" || exit 1

# 끝난 검체 건너뛰기 - 피처 5종이 다 있고 개수가 같으면 끝 (다시 돌리려면 VALREDO=1)
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

# val04_1c.py (JSD 피처) 없으면 복사
if [ ! -f val04_1c.py ]; then
  cp -f "$TRAIN_SRC/06_jsdfeat.py" val04_1c.py 2>/dev/null \
    || echo "!! val04_1c 원본 없음" >> "$L"
fi
# 스크립트 최신화: 검체 폴더의 val*.py 를 학습 원본과 md5 로 맞춤 (config.py 는 안 건드림)
# 학습 스크립트 원본 (config.conf 의 TRAIN_SRC)
_SRCDIR=$(ls -d "$TRAIN_SRC"/*/$v 2>/dev/null | head -1)
[ -n "$_SRCDIR" ] || _SRCDIR=$(ls -d "$TRAIN_SRC" 2>/dev/null | head -1)
if [ -n "$_SRCDIR" ]; then
  # 이름 짝은 config.conf 의 STEP_MAP, 비어 있으면 아래 기본값
  _pairs="${STEP_MAP:-val03_1:01_Preprocessing_for_ML_RawData_Mixing.py
val03_2:02_Preprocessing_for_ML_RawData_Mixing.py
val03_3:03_Preprocessing_for_ML_RawData_Mixing.py
val04_0:04_coverage_check.py
val04_1:05_Feature_Matrix.py
val04_1c:06_jsdfeat.py
val04_1d:07_llrfeat.py}"
  # heredoc 말고 파이프로 읽음
  _MISS=$(mktemp)
  printf '%s\n' "$_pairs" | tr ';' '\n' | while IFS=: read -r _n _f; do
    [ -n "$_n" ] || continue
    _src="$_SRCDIR/$_f"; [ -f "$_src" ] || _src="$(dirname "$_SRCDIR")/$_f"
    # 원본 없으면 _MISS 에 적어 두고 루프 뒤에 중단 (파이프 서브셸이라 exit 안 먹음)
    if [ ! -f "$_src" ]; then
      echo "$_f" >> "$_MISS"
      continue
    fi
    _a=$(md5sum "$_src" | cut -d" " -f1)
    _b=$(md5sum "$_n.py" 2>/dev/null | cut -d" " -f1)
    # md5sum 실패, cp 실패도 _MISS 에 적음 (둘 다 빈 md5 면 같다고 보는 것 방지)
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
  # 단일 리드 코호트도 섞기, 추출 코드가 그대로 처리 (단일 리드 패치 안 씀)
else
  # 원본 폴더 없으면 중단 - 옛 스크립트로 돌면 피처 빠진 채 AUC 0.500
  echo "!! 원본 폴더를 못 찾았습니다. TRAIN_SRC=$TRAIN_SRC" >> "$L"
  echo "   학습 스크립트를 최신화할 수 없습니다. 그대로 두면 옛 스크립트가 돌아" >> "$L"
  echo "   그 피처만 빠진 채 완주합니다 (AUC 0.500). config.conf 의 TRAIN_SRC 를 보십시오." >> "$L"
  exit 1
fi
echo "=== 시작 $(date +%H:%M:%S) $s / $v" >> "$L"
# JSD 칸 목록은 학습이 고른 것 그대로 사용 (BEST.joblib 열과 맞춤)
export JSD_COLS=$METH_ROOT/results/panels/${v%%_*}/$v/step04_ML_classifier/jsd_cols_$v.txt
# 없거나 빈 파일(-s)이면 중단
if [ ! -s "$JSD_COLS" ]; then
  echo "!! 칸 목록이 없습니다. $JSD_COLS" >> "$L"
  echo "   학습(2_train/run_train.sh 의 06)이 이 파일을 만듭니다. 학습을 먼저 돌리십시오." >> "$L"
  echo "   없이 돌리면 검증 자료로 칸을 다시 골라 AUC 0.500 이 납니다." >> "$L"
  exit 1
fi
# val04_1d.py (LLR 피처) 없으면 복사
if [ ! -f val04_1d.py ]; then
  cp -f "$TRAIN_SRC/07_llrfeat.py" val04_1d.py 2>/dev/null     || echo "!! val04_1d 원본 없음" >> "$L"
fi
# LLR 참조표, 후보 풀, JSD 후보, JSD 개수표 - 학습과 같은 GEN 세대 사용
export LLR_REF=${LLR_REF:-$METH_ROOT/results/dmr/j_llrref${GEN}/llrref_K${PATTERN_K:-3}_minr${MIN_READS_PER_WINDOW:-6}.parquet}
export LLR_CAND=${LLR_CAND:-$METH_ROOT/results/dmr/j_candidates${GEN}}
# val04_1c(06_jsdfeat)는 인자 없이 환경변수로만 받음
export JSD_CAND=${JSD_CAND:-$METH_ROOT/results/dmr/j_candidates${GEN}}
export JSD_JCOUNT=${JSD_JCOUNT:-$METH_ROOT/results/dmr/j_jsdcount${GEN}}
[ -f "$LLR_REF" ] || { echo "!! 참조표 없음 $LLR_REF" >> "$L"; exit 1; }
# llr 블록 집합은 학습 mean 표에서 가져옴 (07_llrfeat.py 가 읽음)
export LLR_KEEP_FROM=$METH_ROOT/results/panels/${v%%_*}/$v/step04_ML_classifier
# 검증 실행 표시 - 05_Feature_Matrix 가 블록 버리기를 안 함
export METH_VALIDATE=1
# 참조표는 상위 폴더까지 찍음 (세대마다 파일 이름이 같음)
echo "--- 참조 $(basename $(dirname $LLR_REF))/$(basename $LLR_REF) · 후보 $(basename $LLR_CAND)" >> "$L"
# 도는 순서 = $_pairs (STEP_MAP 또는 기본값) 왼쪽 이름 순서
_STEPS=$(printf '%s\n' "$_pairs" | tr ';' '\n' | cut -d: -f1 | grep -v '^$' | tr '\n' ' ')
[ -n "$_STEPS" ] || { echo "!! 돌 단계 목록이 비었다 — STEP_MAP 을 보십시오" >> "$L"; exit 1; }
echo "--- 단계 $_STEPS" >> "$L"
for st in $_STEPS; do
  # 단계 파일 없으면 중단 (복사 실패)
  [ -f "$st.py" ] || { echo "!! $st.py 가 검체 폴더에 없다 — 복사가 안 됐다" >> "$L"; exit 1; }
  # 메모리 문은 원본 파일명(02_*)으로 판단 (왼쪽 이름은 STEP_MAP 에서 바뀔 수 있음)
  _SRCF=$(printf %s "$_pairs" | tr ';' '\n' | awk -F: -v k="$st" '$1==k{print $2}' | head -1)
  case "$_SRCF" in 02_*) _GATE=1 ;; *) _GATE=0 ;; esac
  # 메모리 문: 섞기 단계만 메모리를 많이 써서 여기만 동시 수 제한
  if [ "$_GATE" = 1 ]; then
    # 잠금 파일 N개(MEMGATE_N)를 토큰으로 씀 - 하나 잡아야 들어감, 죽으면 자동 반납
    #   여유 메모리 MEMGATE_FREE(G) 이상일 때만 시도. 스왑 0 이라 토큰 수 올리지 말 것
    _got=0
    for _i in $(seq 1 540); do
      _a=$(free -g | awk '/^Mem:/{print $7}')
      case "$_a" in ''|*[!0-9]*) _a=0 ;; esac      # 값 이상하면 부족으로 봄
      if [ "$_a" -ge "${MEMGATE_FREE:-22}" ]; then
        for _k in $(seq 1 ${MEMGATE_N:-22}); do
          # 문은 fd 6 사용 (fd 9 는 검체 잠금, 재할당하면 잠금 풀림)
          exec 6>>"$TMP_ROOT/lock/val03_2.$_k"
          flock -n 6 && { _got=1; break; }
        done
      fi
      [ "$_got" = 1 ] && break
      [ "$_i" = 1 ] && echo "--- val03_2 문 대기 (여유 ${_a}G)" >> "$L"
      sleep $((15 + RANDOM % 11))                  # 대기자들 깨는 시각 흩뜨림
    done
    if [ "$_got" != 1 ]; then
      # 3시간 못 잡으면 실패로 나감 (통과시키지 않음, 다음 실행에서 다시 잡힘)
      echo "!! val03_2 문 3시간 초과: 이 검체를 건너뛴다 (나중에 다시 잡힌다)" >> "$L"
      exit 1
    fi
  fi
  echo "--- $st $(date +%H:%M:%S)" >> "$L"
  # fd 6, 9 는 자식에게 안 넘김 (고아 자식이 잠금 쥐고 남는 것 방지)
  python -u "$st.py" 6>&- 9>&- >> "$L" 2>&1 || { echo "!! $st 실패" >> "$L"; exit 1; }
  # 단계 끝나면 문 토큰 바로 반납
  [ "$_GATE" = 1 ] && exec 6>&-         # 문 토큰만 반납, 검체 잠금(fd 9)은 유지
done
echo "=== 끝 $(date +%H:%M:%S)" >> "$L"

# 중간물 정리: 피처 5종이 다 있고 개수 같을 때만 sample_data/panels/<앞머리>/<판>/ 지움
#   results/ 는 안 건드림. 끄려면 VALKEEP=1
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
  # BAM 회수: 디스크 여유가 DISKMIN(G) 아래일 때만 step03 sampled_reads_*.bam 지움
  #   채점은 CSV 만 읽음. 지우면 피처 다시 만들 때 bismark 부터 다시
  #   background_map_*.csv, txt.tar.gz 는 남김
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
