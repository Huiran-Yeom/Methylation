# Methylation-based GBM classifier: panel selection to external validation

Code for building a methylation-pattern classifier for glioblastoma (GBM) from
cfDNA and validating it on an external cohort. The same code, with different
inputs, was run on small-cell lung cancer (SCLC) cell lines as a supporting study.

The tables in `../results/` were produced by commit `e8d5fe5d` (tag
`results-e8d5fe5d`). Later commits change four computations (see
[Reproducibility](#reproducibility)), so check out that tag to regenerate them.

---

## Pipeline

Edit one file, then run four entry points in order.

```bash
# 0) point the config at your data — data roots (METH_ROOT, DATASET_ROOT)
#    and work roots (WORK_ROOT, TMP_ROOT)
vi 0_setup/config.conf

# 1-4)
bash 1_dmr/run_dmr.sh
bash 2_train/run_train.sh
COHORT=<your_cohort> bash 3_validate/run_validate.sh j15bl200_5k_uni 8
COHORT=<your_cohort> bash 4_report/make_report.sh   # same COHORT as stage 3
```

| | Stage | Entry point | What it does |
|---|---|---|---|
| 0 | Setup | `0_setup/config.conf` | paths, depth, ratios, generation |
| 1 | DMR panel selection | `1_dmr/run_dmr.sh` | candidate blocks → JSD → 200-block panel + 6 null panels + pooled BAMs |
| 2 | Training | `2_train/run_train.sh` | in-silico mixing → five features → four models, cross-validated |
| 3 | External validation | `3_validate/run_validate.sh` | mix the external cohort under the training conditions, then score |
| 4 | Reporting | `4_report/make_report.sh` | discrimination tables, paired bootstrap intervals |

Each stage folder has one entry point; its step code sits in `steps/`, numbered
in dependency order. `0_setup/IO.md` lists every step with its inputs and outputs.

### Running every panel and depth

Stage 1 builds all nine panels (3 real + 6 null) in one run and does not depend
on depth. Stages 2 and 3 run one panel at one depth, chosen by `VERSION`, `PANEL`
and `TARGET_DEPTH_READS`. The depth appears both in `VERSION` (`_5k_`) and in
`TARGET_DEPTH_READS`; `0_setup/meth_config.py` stops if they disagree, and
`make_val_config.py` stops if `VERSION` and the panel passed to stage 3 disagree.
Environment variables win over `config.conf`, so this loop edits no file:

```bash
REAL="bl200 jsd200 jsdb200"
NULL="rand1200 rand2200 rand3200 randb1200 randb2200 randb3200"

for d in 5k:5000 50k:50000 100k:100000; do
  tag=${d%%:*}; reads=${d##*:}
  [ "$tag" = 5k ] && panels="$REAL $NULL" || panels="$REAL"   # null panels at 5k only
  for p in $panels; do
    V=j15${p}_${tag}_uni
    export VERSION=$V PANEL=$p TARGET_DEPTH_READS=$reads
    bash 2_train/run_train.sh || continue
    COHORT=<your_cohort> bash 3_validate/run_validate.sh "$V" 8
  done
done
COHORT=<your_cohort> bash 4_report/make_report.sh
```

With the default `LLR_HELDOUT=1` the normal pool is split between the LLR
reference and the training pool, so a sample at 100k usually cannot be filled.
Those cells stop in stage 2 with `풀 부족으로 중단` and the loop moves on; the
reported tables contain 5k and 50k for GBM, plus `jsdb200` at 100k for lung.

On a 32-core server, stage 1 from empty output roots took 59 minutes for GBM and 99 minutes for lung; one stage-2 panel-depth takes about 7 hours; one stage-3 sample
takes about 34 minutes, 33 of them in `bismark`. Stage 3 runs samples in parallel,
up to the memory gate described below.

### Patient plasma cohorts: score without dilution (exploratory)

Stage 3 mixes each validation sample into a normal background at the training
ratios. That is right for tissue or cell-line DNA, which is close to 100% tumour.
Patient plasma is already diluted in the body, and mixing it again at 0.1–5%
pushes the tumour reads below the detection limit. For plasma cohorts, score each
person on their own reads and score same-platform healthy donors the same way as
negatives:

```bash
cp 3_validate/steps/cohorts/patient_plasma.conf.example ~/cohorts/<plasma>.conf   # RATIOS=0.0, 1.0 · NEG_SAMPLES
COHORT=<plasma> bash 3_validate/run_validate.sh j16jsd200_5k_uni 8                  # once per panel
COHORT=<plasma> python 4_report/steps/plasma_baseline.py --cov '<dir>/{sample}.bismark.cov.gz' -o baseline.csv
COHORT=<plasma> python 4_report/steps/plasma_metric.py --baseline baseline.csv j16bl200 j16jsd200 j16jsdb200
```

`plasma_metric.py` reports person-level AUC with a sham control, a label
permutation p-value, and two panel-free baselines (total reads, global mean
methylation) measured on the same permutations. If a baseline alone separates
patients from donors, the result is reported as confounded whatever the panel
AUC. `make_report.sh` refuses such a cohort, because its tables compare samples
against the 0% background rows and mean nothing when donors sit in the sample
slot. Only 5k is practical: one person rarely has 50k reads inside a 200-block
panel. This mode was added after the pre-registered analysis and does not replace
the primary metric (`4_report/steps/primary_metric.py`).

---

## Method in brief

Panel. Among 100 bp blocks covered by several normals, blocks with ≥3 CpGs
become candidates. Within those we count methylation patterns over K=3 CpG
windows and take the 200 blocks with the largest Jensen-Shannon divergence
between cancer and normal (the JSD panel). For comparison we also build a
Baseline panel by moderated t-test plus bootstrap stability ranking.

Six null panels measure how much would be found if blocks were chosen without
regard to cancer.

| | How chosen | Yardstick for |
|---|---|---|
| Null A `rand1-3` | uniform draw from blocks where a K=3 window stands | JSD panel |
| Null B `randb1-3` | stratified draw from the Baseline candidate pool, matching its CpG-count distribution | Baseline panel |

Seeds are fixed at `20260914+i` and `20260914+100+i`.

Five features

| Feature | Definition |
|---|---|
| `mean` | block mean methylation |
| `entropy` | entropy of the within-block beta distribution (10-bin histogram) |
| `jsd` | JSD between the window pattern distribution and the normal reference |
| `pdr` | fraction of fragments whose three CpGs are not all identical (Landau 2014) |
| `llr` | log-likelihood ratio per θ, and its maximum (`LLR_max`) |

`LLR(θ=0) ≡ 0`, so `LLR_max` is floored at zero: a median of 0.00 means "at the
floor", not "missing".

Mixing. Cancer fragments are spiked into a normal background at `MUT_RATIOS`
0, 0.1, 0.5, 1, 2, 2.5, 3, 5, 10 and 100%. The RNG is seeded from `SEED`, the
ratio, the replicate and the sample name, so each combination reproduces and no
two samples share a negative draw. A ratio whose pool cannot be filled is skipped
rather than filled by sampling with replacement, which would make replicates
non-independent.

Feature combinations. Models are cross-validated over all 31 non-empty subsets
of the five features. Two further features (`readent`, `mhl`) were dropped before
results were read because their fill rate at the primary depth was under 50%.

---

## Input data

Not included, for size and access-control reasons. The cohort definition is not
included either: it lists sample identifiers, and for controlled-access datasets
the list may fall under the data access agreement. Structure-only examples:

```
3_validate/steps/cohorts/sample_cohort.conf.example    tissue / cell-line cohorts (mixed)
3_validate/steps/cohorts/patient_plasma.conf.example   patient plasma (not mixed)
  -> copy to ~/cohorts/<cohort>.conf, fill in, add the name to VAL_COHORTS in 0_setup/config.conf
```

| What | Where | Note |
|---|---|---|
| Training normals (BAM) | `NORMAL_BAM_DIR` | public cfDNA |
| Training cancer (cov, BAM) | `GBM_COV_DIR`, `GBM_BAM_DIR` | GBM cell lines |
| Validation cancer BAM | `CANCER` in the cohort definition | EGAD00001003427 (controlled access) |
| Validation normal background | `NORMAL` / `NORMALS_*` in the same file | zero overlap with training |

The GBM validation samples are tissue DNA, not plasma cfDNA. This simulates
clinical validation; confounding by the tissue/cfDNA difference cannot be
excluded.

### What a fresh clone still needs

| Needed | Where the code looks | If missing |
|---|---|---|
| Normal cov files | `NORMAL_SET`, default `$METH_ROOT/sample_data/_input/normal_pub15/*.cov.gz` — under the output root, not `DATASET_ROOT` | step 01: `정상 cov 가 없다` |
| Bisulfite genome index | `$METH_ROOT/data/Bisulfite_Genome/` | step 03 runs bismark without `--genome_folder`, which changes its behaviour silently |

Step `10c` builds the three panel-union CSVs and the two sample lists that step
11 consumes.

Data folders. Changing `METH_ROOT` and `DATASET_ROOT` moves every path. If your
layout under those roots differs, set these instead of editing the steps (each
value is read in two or three places):

| Variable | Default | Read by |
|---|---|---|
| `NORMAL_SET` | `$METH_ROOT/sample_data/_input/normal_pub15` | steps 01, 02, 03 |
| `GBM_COV_DIR` | `$DATASET_ROOT/training/GBM_cell-line/cov` | steps 01, 02, 03 |
| `NORMAL_BAM_DIR` | `$DATASET_ROOT/for_in_silico_test/normal_cfDNA_public/aligned_bam` | steps 03, 11 |
| `GBM_BAM_DIR` | `$DATASET_ROOT/training/GBM_cell-line/aligned_bam` | steps 03, 11 |

BAM and cov files are paired by the sample id. Step 03 checks every pair before
it starts and stops if any is missing:

```
  GBM      짝 30개 · 못 찾음 0개
  Normal   짝 0개 · 못 찾음 15개
중단: 짝 BAM 을 못 찾은 검체가 15개 있습니다.
```

Any missing pair is fatal, because step 01 writes `samples.csv` from the cov list
while step 03 scans only what it could pair, and step 04 derives its coverage
threshold from `samples.csv`. `ALLOW_SKIP=1` proceeds anyway and says so.
`--sample <id>` runs one sample.

Naming rules. These variables adapt sample ids and grouping to another dataset
without renaming files. Their defaults reproduce the GBM results.

| Variable | Default | What it does |
|---|---|---|
| `SID_RE` | first `_`-separated field | sample id from the filename (first capture group) |
| `GROUP_RE` | `SNU<n>-<rep>` → `SNU<n>` | groups replicate files into one independent unit, cancer side only |
| `PERSON_RE` | `pat_1_…` / `ind_3_…` | stage 4 only: groups validation samples by person |
| `DEPTHS` | `5k 15k 50k 100k` | panel depths read in stages 3 and 4 |
| `CANCER_LABEL` | `GBM` | the `type` value stage 1 writes in `samples.csv` and later steps compare against |

Steps 01 and 02 stop if two files resolve to the same sample id. `GROUP_RE` is
not applied to normals.

### The lung (SCLC) run

The lung tables used the same code with `GEN=16`, `CANCER_LABEL=SCLC`, the cell
line cov/BAM folders in `GBM_COV_DIR` / `GBM_BAM_DIR`, the training normals in
`NORMAL_SET` / `NORMAL_BAM_DIR`, and

```bash
SID_RE='(.+?)__' GROUP_RE='(.+?)_\d+d$'
```

so that `H1876_0d__SRR…` and `H1876__SRR…` (one cell line sequenced twice) become
ids `H1876_0d` and `H1876` grouped into one unit.

Stages 2 and 4 take the normal/cancer split from the `type` column of
`samples.csv`. Inside `2_train`, sample ids are built as `Normal1…` / `GBM1…`
whatever `CANCER_LABEL` is, and the classifier keys on that prefix.

Pooled BAMs. A pooled BAM holds only reads inside its own panel regions, so step
11 builds three sets (`_cellline`, `_rand`, `_randb`) and `2_train/steps/config.py`
picks the one matching `PANEL`.

Memory. Stage 3 waits before the mixing step until `free -g` reports at least
22 GB available (about 1.35 GB per sample) and gives up after about 3 hours. At
most `MEMGATE_N` samples (default 22) hold a token at once. Both are environment
variables: `MEMGATE_FREE=<GB>`, `MEMGATE_N=<n>`.

---

## Environment

Versions used are pinned in `requirements.txt`; `0_setup/ENVIRONMENT.md` explains
why each matters. Linux only: the pipeline uses `flock`, `/proc`, `bismark`,
`bowtie2` and `samtools`.

```bash
conda create -n meth python=3.9
conda activate meth
pip install -r requirements.txt
# bismark, bowtie2, samtools are installed separately: see ENVIRONMENT.md
```

---

## Reproducibility

Which commit made which numbers:

| Numbers | Commit |
|---|---|
| `../results/` (GBM and lung, training CV and validation) | `e8d5fe5d`, tag `results-e8d5fe5d` |
| this branch | `e8d5fe5d` plus four changes from `79979d24`: the held-out split in `05_llrref` uses the whole `samples.csv`; validation llr blocks come from the training mean table (`LLR_KEEP_FROM`); validation mode is also signalled by `METH_VALIDATE=1`; in stage 2 the threshold predictions are imputed inside each fold and Logistic Regression uses `max_iter=2000`. These change `cv_all.csv` and the thresholds, so this branch does not regenerate `../results/`. |

A from-scratch re-run of `e8d5fe5d` (all of stage 1 on empty output roots, then
stages 2–3, BLAS threads pinned to 1) was compared with the stored outputs at a
tolerance fixed before the run (1e-9 on CV tables and per-row scores). It did
not pass that tolerance:

- stage 1 matched: nine panels and the LLR reference byte for byte, pooled BAMs
  by read count;
- in stage 2, Random Forest, SVM and XGBoost CV values matched exactly; 4–10% of
  the Logistic Regression rows differed (at most 0.0135), and two of 403 lung
  BEST models switched between Logistic Regression and XGBoost;
- the primary metric matched to four decimals in all nine cells.

The cause is the BLAS thread count. On one cell, re-running step 08 on the same
feature CSVs gave bit-identical results run to run; with threads pinned to 1 it
matched the re-run, and with threads unpinned it matched the stored tables in all
1,116 rows. The early run scripts that trained the stored tables did not pin
threads (checked cell by cell for this one cell only). Logistic
Regression on `main` (`max_iter=2000`, no convergence warnings) still differs
between the two thread settings (at most 0.0005), so record the thread setting
as part of any reproduction. To regenerate `../results/`, leave
`OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS` and `MKL_NUM_THREADS` unset.

`../results/main/` holds the same tables re-run on `main`; the primary metric
moves by at most 0.002 and every cell stays below 0.70.

An earlier check, on the code as of 2026-09-30 (before the held-out LLR change),
re-ran parts of the pipeline and matched the stored outputs byte for byte:

| Stage | Re-run | Result |
|---|---|---|
| 1 | all nine panels, from candidate blocks | `DMR_confirmed_*.csv` md5 9 / 9 |
| 2 | `bl200`, `jsd200`, `jsdb200` at 5k | `BEST.joblib` md5 31 / 31 per panel |
| 3 | `bl200` at 5k, 178 samples | feature CSVs 6,230 / 6,230, `y_prob_all_*.csv` 178 / 178 |
| 3 | `jsd200`, `jsdb200` at 5k, 4 samples each | feature CSVs 140 / 140, `y_prob_all_*.csv` 4 / 4 |

A fixed seed is not a guarantee. Numbers can move with library versions
(especially `scikit-learn`), BLAS thread count and `n_jobs`, `bismark` /
`bowtie2` / `samtools` versions, and the OS. `JSD_LIMIT` caps how many reads
stage 1 scans; use it to check that a stage runs, not to compare values.

Two stages skip work that already exists: mixing skips finished ratios and
scoring skips samples whose model fingerprint matches. Re-run into empty roots
when comparing.

---

## Things to know when reading the code

Config lives in one place. Edit `0_setup/config.conf`; environment variables win
over it. Each output folder gets a `config_snapshot.txt` recording where every
value came from (env / conf / default).

File names are coupled to the config. Validation copies the training scripts
into each sample's folder; if you rename anything under `2_train/steps/`, update
`STEP_MAP` in `0_setup/config.conf`. Otherwise that feature is missing, scoring
fills it with the training mean, and you get AUC 0.500. A guard in the copy loop
stops this.

`GEN` is the suffix on panel output folders (15 for GBM, 16 for lung). Several
guards stop if training and validation see different generations.

Scripts refuse to be imported, because importing an unguarded script runs its
body. Run them as `python <file>`.

Stage 1 refuses to write into a non-empty output folder, except for steps that
add to a folder an earlier step made (02, 04b, 04c, 09, 10). To overwrite, retype
the folder name: `OVERWRITE=<folder>`. Stages 2–4 have no such guard.

Locks and logs live under `TMP_ROOT` (default `~/tmp`).

Output names are ASCII. `0_setup/RENAMES.txt` maps the script names used in the
paper's methods to the names here. Code comments are in Korean.

---

## Layout

```
README.md  requirements.txt  .gitignore  .gitattributes
0_setup/      config.conf (the file you edit)  meth_config.py  check_fixes.py
              IO.md  ENVIRONMENT.md  RENAMES.txt
1_dmr/        run_dmr.sh   steps/ 01_candidates 02_blockcpg 03_jsdcount 04_jsd
                                  04b_jsdnull 04c_jsdboot 05_llrref 06_select
                                  06b_select_supplement 07_panel 08_export_panels
                                  09_null_a 10_null_b 10c_panel_union 11_pool
                                  run_pool.sh moderated_t
2_train/      run_train.sh steps/ 01-03_Preprocessing…Mixing 04_coverage_check
                                  05_Feature_Matrix 06_jsdfeat 07_llrfeat
                                  08_ML_classifier 09_save_best_j config
3_validate/   run_validate.sh steps/ run_one_sample make_val_config score_one_sample
                                  paths val_single
                                  cohorts/ sample_cohort.conf.example patient_plasma.conf.example
4_report/     make_report.sh steps/ make_tables primary_metric plasma_metric plasma_baseline
5_analysis/   one-off analyses behind the pre-registration addenda (see its README)
```
