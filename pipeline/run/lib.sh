#!/bin/bash
# 드라이버 공용 함수(2026-10-27 정리판). 옛 jobs_lib.sh · r1_env.sh · r2_env.sh 를 한 파일로 합쳤다.
#   옛 jobs_lib.sh 와 다른 점:
#     - 코드 경로 P 를 정하지 않는다(드라이버가 정한다 · 옛 판은 P 를 Methylation-fix 로 덮어 불러오는 순서에 따라 옛 코드로 돌았다).
#     - train_cell 의 관문을 세게 했다: mut_* 폴더가 하나도 없으면 실패, BEST 모델은 조합 수(EXPECT_BEST=31 고정)를 다 채워야 성공.
#   env_luad2 가 env_luad 를 부르므로 둘을 같은 파일에 둔다(옛 판은 r1_env 를 먼저 불러야 했다).
export PATH=$HOME/anaconda3/envs/myenv/bin:$PATH
export BISMARK_MULTICORE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
NSLICE=${NSLICE:-4}
EXPECT_BEST=31        # 5특징의 비공집합 31조합 — 밖에서 덮지 못하게 고정한다

env_lung(){   # env_lung <METH_ROOT> <작업 이름> [cand]
  unset USE_NORMAL USE_GBM JSD_CAND LLR_CAND
  export GEN=16 CANCER_LABEL=SCLC NSMP=8
  export NORMAL_SET=/ssd_data/dataset/Normal/Train/cov NORMAL_BAM_DIR=/ssd_data/dataset/Normal/Train/bam
  export GBM_COV_DIR=/ssd_data/dataset_for_lung/stage1_catalog/tumor_cellline/cov
  export GBM_BAM_DIR=/ssd_data/dataset_for_lung/stage1_catalog/tumor_cellline/aligned_bam
  export SID_RE='(.+?)__' GROUP_RE='(.+?)_\d+d$'
  export METH_ROOT=$1 WORK_ROOT=$HOME/work61/$2 TMP_ROOT=$HOME/tmp/t61/$2
  if [ "${3:-}" = cand ]; then
    export JSD_CAND=/ssd_data/Methylation_fix/lung/results/dmr/j_candidates16
    export LLR_CAND=/ssd_data/Methylation_fix/lung/results/dmr/j_candidates16
  fi
  mkdir -p "$WORK_ROOT" "$TMP_ROOT"
}
env_gbm(){    # env_gbm <METH_ROOT> <작업 이름>
  unset NORMAL_SET GBM_COV_DIR GBM_BAM_DIR SID_RE GROUP_RE JSD_CAND LLR_CAND
  export GEN=15 CANCER_LABEL=GBM NSMP=12
  export USE_NORMAL=$HOME/tmp/use15_normal.txt USE_GBM=$HOME/tmp/use15_gbm.txt
  export NORMAL_BAM_DIR=/ssd_data/dataset/Normal/Train/bam
  export METH_ROOT=$1 WORK_ROOT=$HOME/work61/$2 TMP_ROOT=$HOME/tmp/t61/$2
  mkdir -p "$WORK_ROOT" "$TMP_ROOT"
}

# ── 다음 판 1안 · 2안 (LUAD) ─────────────────────────────────────────────
R1_RAW=/ssd_data/dataset_for_lung/stage1_catalog/tumor_cellline_luad
R2_TISSUE=/ssd_data/dataset_for_lung/stage1_catalog/normal_lung_tissue
r1_root(){ case "$1" in 10) echo /ssd_data/Methylation_luad_a1;; 5) echo /ssd_data/Methylation_luad_a1u;; *) return 1;; esac; }
r2_root(){ case "$1:$2" in 10:t) echo /ssd_data/Methylation_luad_a2t;; 5:t) echo /ssd_data/Methylation_luad_a2tu;;
                           10:m) echo /ssd_data/Methylation_luad_a2;;  5:m) echo /ssd_data/Methylation_luad_a2u;; *) return 1;; esac; }
env_luad(){   # env_luad <판 10|5>
  local s=$1 R; R=$(r1_root "$s") || { echo "판은 10 또는 5"; return 1; }
  unset USE_NORMAL USE_GBM JSD_CAND LLR_CAND COLSET_FROM SCORE_UNIT CONVERGENCE_GATE
  export GEN=17 CANCER_LABEL=LUAD NSMP=8
  export NORMAL_SET=/ssd_data/dataset/Normal/Train/cov NORMAL_BAM_DIR=/ssd_data/dataset/Normal/Train/bam
  if [ "$s" = 10 ]; then export GBM_COV_DIR=$R1_RAW/cov GBM_BAM_DIR=$R1_RAW/aligned_bam
  else export GBM_COV_DIR=$R/inputs/cov GBM_BAM_DIR=$R/inputs/bam; fi
  export SID_RE='(.+?)_1_val_1_' GROUP_RE='(A549|H358)_'
  export LLR_HELDOUT=0 BLAS_THREADS=1
  export PANELS=bl200,jsd200
  export METH_ROOT=$R WORK_ROOT=$HOME/work61/r1_luad$s TMP_ROOT=$HOME/tmp/t61/r1_luad$s
  mkdir -p "$WORK_ROOT" "$TMP_ROOT"
}
env_luad2(){   # env_luad2 <판 10|5> <t|m>
  local s=$1 role=$2 R; R=$(r2_root "$s" "$role") || { echo "판 10|5 · 역할 t|m"; return 1; }
  env_luad "$s" || return 1
  [ "$role" = t ] && export NORMAL_SET=$R2_TISSUE/cov NORMAL_BAM_DIR=$R2_TISSUE/aligned_bam
  export METH_ROOT=$R WORK_ROOT=$HOME/work61/r2_luad${s}$role TMP_ROOT=$HOME/tmp/t61/r2_luad${s}$role
  mkdir -p "$WORK_ROOT" "$TMP_ROOT"
}

best_n(){      # best_n <판 버전 폴더> — BEST 모델 수
  ls "$1"/step04_ML_classifier/models/*/pooled/BEST.joblib 2>/dev/null | wc -l
}

train_cell(){  # train_cell <판> <깊이> <리드 수> — 드라이버가 P · env_* 를 먼저 정한다
  local p=$1 d=$2 n=$3 rc=0 i s T0 dir a b nd=0
  : "${P:?P(코드 경로)를 먼저 정한다}"
  export PANEL=$p VERSION=j${GEN}${p}_${d}_uni TARGET_DEPTH_READS=$n
  local l=$TMP_ROOT/train_${p}_${d}.log pv=j${GEN}${p}
  echo "[$(date '+%m-%d %H:%M:%S')] $VERSION · METH_ROOT=$METH_ROOT · 03 조각 $NSLICE · 코드 $P" > "$l"
  T0=$(date +%s)
  for s in 01 02; do [ "$rc" = 0 ] && { bash "$P/2_train/run_train.sh" $s >> "$l" 2>&1 || rc=1; }; done
  if [ "$rc" = 0 ]; then
    local pids=()   # wait -n 은 같은 셸의 다른 배경 작업도 센다 — 조각 PID 로 하나씩
    for i in $(seq 1 "$NSLICE"); do SLICE_N=$NSLICE SLICE_I=$i bash "$P/2_train/run_train.sh" 03 > "$l.s$i" 2>&1 & pids+=($!); done
    for i in "${pids[@]}"; do wait "$i" || rc=1; done
    for i in $(seq 1 "$NSLICE"); do echo "── 03 조각 $i/$NSLICE" >> "$l"; cat "$l.s$i" >> "$l"; rm -f "$l.s$i"; done
    for dir in "$METH_ROOT/results/panels/$pv/$VERSION/step03_Preprocessing_for_ML"/mut_*_reads; do
      [ -d "$dir" ] || continue; nd=$((nd + 1))
      a=$(ls "$dir"/sampled_reads_*.bam 2>/dev/null | wc -l)
      b=$(ls "$METH_ROOT/sample_data/panels/$pv/$VERSION/step04_ML_classifier/$(basename "$dir")"/sampled_reads_*.bismark.cov.gz 2>/dev/null | wc -l)
      [ "$a" -gt 0 ] && [ "$a" -eq "$b" ] || { echo "!! 조각 관문: $(basename "$dir") BAM $a · cov $b" >> "$l"; rc=1; }
    done
    [ "$nd" -gt 0 ] || { echo "!! 조각 관문: mut_*_reads 폴더가 없다" >> "$l"; rc=1; }
  fi
  for s in 04 05 06 07 08 09; do [ "$rc" = 0 ] && { bash "$P/2_train/run_train.sh" $s >> "$l" 2>&1 || rc=1; }; done
  if [ "$rc" = 0 ]; then
    i=$(best_n "$METH_ROOT/results/panels/$pv/$VERSION")
    [ "$i" = "$EXPECT_BEST" ] || { echo "!! BEST 모델 $i 개(${EXPECT_BEST} 여야 함)" >> "$l"; rc=1; }
  fi
  echo "[$(date '+%m-%d %H:%M:%S')] 끝 rc=$rc · $(( ($(date +%s)-T0)/60 ))분" >> "$l"
  return $rc
}
