# 단계별 입력·출력

`<뿌리>` 는 `config.conf` 의 `METH_ROOT`, `<자료>` 는 `DATASET_ROOT` 입니다.
`<판>` 은 `VERSION`(예 `j15bl200_5k_uni`), `<패널>` 은 `PANEL`(예 `bl200`)입니다.

---

## 1단계 · 패널 선정 → `1_dmr/steps/`

번호가 곧 실행 순서입니다. `<G>` 는 `config.conf` 의 `GEN`(지금 `15`)입니다.

| 스크립트 | 입력 | 출력 |
|---|---|---|
| `01_candidates.py` | 정상 cov + `<자료>/training/GBM_cell-line/cov` | `<뿌리>/results/dmr/j_candidates<G>` |
| `02_blockcpg.py` | 위 후보 | 같은 폴더의 `blocks_cpg` |
| `03_jsdcount.py` | 후보 + `<자료>/…/aligned_bam` | `j_jsdcount<G>/count_*_K<PATTERN_K>.parquet` |
| `04_jsd.py` | `j_jsdcount<G>` | `j_jsd<G>/jsd_K<K>_minr<MINR>.parquet` |
| `04b_jsdnull.py` | `j_jsdcount<G>` | `j_jsd<G>/null_K<K>_minr<MINR>.parquet` — 순열 귀무 분포 |
| `04c_jsdboot.py` | `j_jsdcount<G>` | `j_jsd<G>/boot_K<K>_minr<MINR>.parquet` — 복원추출 100회 순위 빈도 |
| `05_llrref.py` | `j_jsdcount<G>` + 후보 | `j_llrref<G>/llrref_K<K>_minr<MINR>.parquet` |
| `06_select.py` | 후보 + `j_jsd<G>` | `j_panel<G>` — JSD 상위 200 (jsd200) |
| `06b_select_supplement.py` | `boot_*` + `null_*` + 후보 | `j_panel<G>_supplement/panel_jsd_supplement_n200.csv` — **jsdb200** |
| `07_panel.py` | 후보 | `j_panel_dmr<G>/bl200_d10_cellline` — Baseline 200 (bl200) |
| `08_export_panels.py` | 06·06b·07 의 셋 | `j_panel_dmr<G>_panels/<패널>/DMR_confirmed_<패널>.csv` |
| `09_null_a.py` | 후보 | 귀무A 3판 (`rand1-3`) · 씨앗 `20260914+i` |
| `10_null_b.py` | Baseline 모집단 | 귀무B 3판 (`randb1-3`) · 씨앗 `20260914+100+i` |
| `11_pool.py` (`run_pool.sh`) | 세포주·정상 BAM + `panel_union3_cellline.csv` | `j_pool<G>_cellline/pool_{gbm,normal}_all.bam` |
| `moderated_t.py` | — | 06·06b·07 이 import 하는 조정 t검정·BH 보정 |

**이 단계의 최종 산출물**은 둘입니다.

1. 9판의 `DMR_confirmed_<패널>.csv` (실제 3 + 귀무 6) — 2단계가 패널로 읽습니다.
2. `j_pool<G>_cellline/pool_{gbm,normal}_all.bam` — 2단계 step 01 이 **섞기의 원천
   풀**로 읽습니다(`config.py` 의 `FULL_GBM_BAM`·`FULL_NORMAL_BAM`).

> 2026-09-29: 04b·04c·06b·11 을 「감도용·미사용」으로 보고 빼려다 되돌렸습니다.
> 04b·04c 는 06b 의 입력이고 06b 는 보고하는 세 판 중 `jsdb200` 을 만듭니다.
> 11 의 풀 BAM 은 2단계 첫 걸음이 읽습니다. 넷 다 본선입니다.
> 다만 `run_pool.sh` 의 입력 `panel_union3_cellline.csv` 를 만드는 코드는 여기
> 없습니다 — 1단계 패널 산출물을 손으로 합친 것입니다(README 에 적어 두었습니다).

---

## 2단계 · 학습 → `2_train/steps/`

공통 입력은 `config.py` 가 읽는 `config.conf` 입니다.
산출물은 전부 `<뿌리>/results/panels/<접두>/<판>/` 아래입니다(아래 `RS`).
중간물은 `<뿌리>/sample_data/panels/<접두>/<판>/` 아래입니다(아래 `SD`).

| 스크립트 | 입력 | 출력 |
|---|---|---|
| `01_…Mixing.py` | 풀 BAM 두 개 + `DMR_confirmed_<패널>.csv` | `SD/step03_…/mut_sorted_<판>.bam`<br>`wt_sorted_<판>.bam` + 이름 목록 |
| `02_…Mixing.py` | 위 정렬 BAM | `RS/step03_…/mut_*_reads/sampled_reads_*.bam`<br>`mixing_plan_<모드>_<판>.csv` |
| `03_…Mixing.py` | 섞은 BAM | bismark 정렬 결과 `txt.tar.gz` |
| `04_coverage_check.py` | 위 | 커버리지 점검 로그 |
| `05_Feature_Matrix.py` | bismark cov | `RS/step04_ML_classifier/<n>_mean_<판>.csv`<br>`<n>_entropy_<판>.csv` |
| `06_jsdfeat.py` | 섞은 BAM + `j_jsdcount` 기준 분포 | `<n>_jsd_<판>.csv` · `<n>_pdr_<판>.csv` |
| `07_llrfeat.py` | 섞은 BAM + `llrref_K<K>_minr<MINR>.parquet` | `<n>_llr_<판>.csv` |
| `08_ML_classifier.py` | 피처 5종 | `cv_<n>_<판>.csv` (조합 31 × 모델 4) |
| `09_save_best_j.py` | `cv_*.csv` + 피처 | `RS/step04_ML_classifier/models/<조합>/pooled/BEST.joblib` |

`<n>` 은 비율의 종양 리드 수입니다(예 `25` = 0.5%).
`BEST.joblib` 에는 모델과 문턱 `thr_spec95` 가 같이 들어갑니다.

---

## 3단계 · 검증 → `3_validate/steps/`

코호트 설정은 `3_validate/steps/cohorts/*.conf` 이고, 여기서 **검증 음성 배경과 반복 수**를 정합니다.

| 스크립트 | 입력 | 출력 |
|---|---|---|
| `run_one_sample.sh <검체> <판>` | 환자 조직 BAM + 코호트가 가리키는 정상 배경 BAM<br>2단계의 `BEST.joblib` | `<뿌리>/<코호트>_val/<검체>/results/<판>/step04_ML_classifier/`<br>피처 5종 CSV |
| `score_one_sample.py` | 피처 5종 + `BEST.joblib` | 같은 폴더의 `y_prob_all_<판>_all_combos.csv` |

**검체마다 나오는 행 수** = 조합 × 농도 × 반복 × (양성·음성).

`run_one_sample.sh` 는 채점을 하지 않습니다 — 피처까지만 만듭니다. 검체 목록을 돌며
`run_one_sample.sh` 를 부르고, 그 다음 검체마다 `score_one_sample.py` 로 채점하는 것은
진입점 `3_validate/run_validate.sh` 가 합니다. 그래서 이 단계에서 직접 부르는 것은
진입점 하나뿐입니다.

> 2026-09-29: 예전에는 여기에 병렬·보충 런처(`run_all_samples.sh` · `run_missing.sh` ·
> `score_parallel.sh`)와 종합표(`aggregate_scores.py`)가 따로 있었습니다. 운영 편의용이고
> 결과에는 영향이 없어 이 저장소에서는 뺐습니다 — 진입점이 그 일을 직접 합니다.
> (표는 검체별 채점 파일을 직접 읽으므로 종합표를 쓰지 않습니다.)

---

## 4단계 · 표 → `4_report/steps/`

| 스크립트 | 입력 | 출력 |
|---|---|---|
| `make_tables.py <판>` | `y_prob_all_*_all_combos.csv` 전부 | 표 텍스트(표준출력) — 조합별 AUC, 검체 단위와 사람 단위 |
| `primary_metric.py` | 같음 | 1차 지표 AUC + 짝지은 부트스트랩 구간 + 특이도 |

`make_tables.py` 는 `COHORT` 환경변수로 코호트를 고릅니다.
`primary_metric.py` 는 `--depth` · `--ratio` · `--combo` · `--nboot` 을 받습니다.

---

## 운영·감시 — 이 저장소에 없습니다

메모리 감시·진행 현황·이어달리기 드라이버는 결과를 만들지 않으므로 넣지 않았습니다.
그 스크립트들이 읽던 자리(로그·잠금)는 `config.conf` 의 `TMP_ROOT` 아래입니다.

## 리드 수를 바꿀 때

`config.conf` 에서 **두 곳**을 같이 고칩니다.

```
VERSION=j15bl200_50k_uni      이름 안의 5k 를 50k 로
TARGET_DEPTH_READS=50000      실제 리드 수
```

이름과 값이 따로 돌기 때문에 한쪽만 고치면 **폴더 이름은 50k 인데 내용은 5k** 가 됩니다.
