# 진입점 공통 머리 (source 해서 씀)
#   1) 이 파일 옆 config.conf 읽기. 환경변수가 우선 (meth_config.py 와 같은 규칙)
#   2) BLAS 스레드 고정
#   3) 암 쪽 자료 이름: 새 이름 CANCER_* 과 옛 이름 GBM_* 둘 다 받음
# 진입점은 0_setup 을 찾아 이 파일을 부름. 트리 밖에서 돌리면 METH_CONF_DIR 로 0_setup 을 알려 줌

# conf 키를 미리 보관해 뒀다가 conf 읽은 뒤 되돌림
_CONFFILE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/config.conf"
_PREKEYS=$(awk -F= '/^[A-Za-z_][A-Za-z0-9_]*=/{print $1}' "$_CONFFILE" 2>/dev/null | tr '\n' ' ')
for _k in $_PREKEYS; do eval "_pre_$_k=\${$_k-}"; done
[ -f "$_CONFFILE" ] || { echo "멈춤: $_CONFFILE 가 없다 (load_config.sh 옆에 config.conf 가 있어야 함)" >&2; exit 1; }
set -a; . "$_CONFFILE"; METH_CONF_DIR=$(dirname "$_CONFFILE"); set +a
for _k in $_PREKEYS; do
  eval "_v=\$_pre_$_k"
  if [ -n "${_v:-}" ]; then eval "export $_k=\$_pre_$_k"; fi
done
unset _v
: "${METH_CONF_DIR:?0_setup/config.conf 를 못 찾았다 — 트리 밖에서 돌리려면 METH_CONF_DIR 를 직접 주라}"

# BLAS 스레드 고정 - 로지스틱 회귀 값이 스레드 수 따라 조금 달라짐
# 환경변수에 이미 있으면 그 값 사용, BLAS_THREADS=none 이면 고정 안 함
# 진입점(run_*.sh, make_report.sh)에서만 고정. steps/ 직접 실행하면 안 걸림
_bt=${BLAS_THREADS:-1}
case "$_bt" in
  [Nn][Oo][Nn][Ee]) ;;
  ''|*[!0-9]*|0) echo "멈춤: BLAS_THREADS 는 none 또는 1 이상의 정수여야 합니다: $_bt" >&2; exit 1 ;;
  *) export OMP_NUM_THREADS=${OMP_NUM_THREADS:-$_bt} OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-$_bt} MKL_NUM_THREADS=${MKL_NUM_THREADS:-$_bt}
     [ "$OMP_NUM_THREADS" = "$OPENBLAS_NUM_THREADS" ] && [ "$OMP_NUM_THREADS" = "$MKL_NUM_THREADS" ]        || echo "주의: 스레드 값이 갈린다 · OMP=$OMP_NUM_THREADS OPENBLAS=$OPENBLAS_NUM_THREADS MKL=$MKL_NUM_THREADS" >&2 ;;
esac

# 암 쪽 자료 이름 맞추기: 하나만 있으면 둘 다 그 값, 둘 다 있고 다르면 멈춤
_same_name(){   # _same_name <새 이름> <옛 이름>
  local a=${!1-} b=${!2-}
  if [ -n "$a" ] && [ -n "$b" ] && [ "$a" != "$b" ]; then
    echo "멈춤: $1 와 $2 가 다릅니다 ($a · $b). 하나만 주십시오" >&2; exit 1
  fi
  if [ -n "$a$b" ]; then export "$1=${a:-$b}" "$2=${a:-$b}"; fi
  return 0
}
_same_name CANCER_COV_DIR GBM_COV_DIR
_same_name CANCER_BAM_DIR GBM_BAM_DIR
_same_name USE_CANCER USE_GBM
