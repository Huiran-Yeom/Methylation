# Methylation-based GBM classifier: panel selection to external validation

Code for building a methylation-pattern classifier for glioblastoma (GBM) from
cfDNA and validating it on an external cohort. This folder contains the code that
produced the numbers reported in the paper, unmerged and unrewritten.

---

## Pipeline: five commands

You edit one file and run four. Everything else is a step invoked for you.

```bash
# 0) point the config at your data (four paths) — the only file you edit
vi 0_setup/config.conf

# 1-4) in order
bash 1_dmr/run_dmr.sh
bash 2_train/run_train.sh
COHORT=<your_cohort> bash 3_validate/run_validate.sh j15bl200_5k_uni 8
COHORT=<your_cohort> bash 4_report/make_report.sh   # same COHORT as stage 3
```

| | Stage | Entry point | What it does |
|---|---|---|---|
| **0** | Setup | `0_setup/config.conf` | Paths, depth, ratios, generation |
| **1** | DMR panel selection | `1_dmr/run_dmr.sh` | candidate blocks → JSD → 200-block panel + 6 null panels + pooled BAMs |
| **2** | Training | `2_train/run_train.sh` | in-silico mixing → five features → four models, cross-validated |
| **3** | External validation | `3_validate/run_validate.sh` | mix the external cohort under identical conditions, then score |
| **4** | Reporting | `4_report/make_report.sh` | discrimination tables · paired bootstrap intervals |

Each stage folder holds exactly **one entry point**; its step code sits in
`steps/`, numbered in dependency order. `0_setup/IO.md` lists every step with its
inputs and outputs; you do not need it to run the pipeline.

### Running all fifteen panel-depth combinations

**Stage 4 reports nine panels at every depth; stages 2 and 3 run one panel at one
depth.** Stage 1 builds all nine panels (3 real + 6 null) in a single run and does
not depend on depth. Training and validation work on whichever `VERSION`/`PANEL`
you point them at, so the full report needs stages 2–3 run once per panel *and*
per depth. The reported tables are **15 combinations**: the three real panels at
5k / 50k / 100k, and the six null panels at 5k only.

Depth lives in two places that must agree: the `5k` inside `VERSION` and
`TARGET_DEPTH_READS`. Both are checked. `0_setup/meth_config.py` stops if they
disagree, and `make_val_config.py` stops if `VERSION` and the panel you pass to
stage 3 disagree. Environment variables win over `config.conf`, so this loop edits
no file:

```bash
REAL="bl200 jsd200 jsdb200"
NULL="rand1200 rand2200 rand3200 randb1200 randb2200 randb3200"

for d in 5k:5000 50k:50000 100k:100000; do
  tag=${d%%:*}; reads=${d##*:}
  # null panels were reported at 5k only
  [ "$tag" = 5k ] && panels="$REAL $NULL" || panels="$REAL"
  for p in $panels; do
    V=j15${p}_${tag}_uni
    export VERSION=$V PANEL=$p TARGET_DEPTH_READS=$reads
    bash 2_train/run_train.sh
    COHORT=<your_cohort> bash 3_validate/run_validate.sh "$V" 8
  done
done
COHORT=<your_cohort> bash 4_report/make_report.sh          # all nine panels, all depths
```

Measured on a 32-core server: stage 1 takes 26 minutes for all nine panels; one
stage-2 panel-depth takes about 7 hours; one stage-3 sample takes about 34
minutes, of which 33 are `bismark`. Stage 3 parallelises across samples, but the
memory gate (below) caps how many can be in the mixing step at once. Panels are
independent, so running several at the same time costs little — bismark spends
most of its wall time waiting rather than computing.

---

## Method in brief

**Panel.** Among 100 bp blocks covered by several normals, blocks with ≥3 CpGs
become candidates. Within those, we count the distribution of methylation
*patterns* over K=3 CpG windows and take the 200 blocks with the largest
Jensen-Shannon divergence between cancer and normal (the JSD panel). For
comparison we also build a Baseline panel by moderated t-test plus bootstrap
stability ranking.

**Six null panels.** These measure "how much would we get if blocks were chosen
*without regard to cancer*."

| | How chosen | Yardstick for |
|---|---|---|
| Null A `rand1-3` | uniform draw from blocks where a K=3 window stands | JSD panel |
| Null B `randb1-3` | stratified draw from the Baseline candidate pool, matching the Baseline CpG-count distribution | Baseline panel |

Seeds are fixed at `20260914+i` and `20260914+100+i`, so both reproduce exactly.

**Five features**

| Feature | Definition |
|---|---|
| `mean` | block mean methylation |
| `entropy` | entropy of the within-block beta distribution (10-bin histogram) |
| `jsd` | JSD between the window pattern distribution and the normal reference |
| `pdr` | fraction of fragments whose three CpGs are not all identical (Landau 2014) |
| `llr` | log-likelihood ratio per θ, and its maximum (`LLR_max`) |

Because `LLR(θ=0) ≡ 0`, `LLR_max` is floored at zero. **A median of 0.00 means
"at the floor", not "missing".**

**Mixing.** Cancer fragments are spiked into a normal background. The shipped
`MUT_RATIOS` are 0, 0.1, 0.5, 1, 2, 2.5, 3, 5, 10 and 100%. The RNG is seeded
from `SEED`, the ratio, the replicate **and the sample name**, so each combination
reproduces and no two samples share a negative draw. Ratios whose pool cannot be
filled are skipped rather than **filled by sampling with replacement**: filling
would make the replicates non-independent and inflate apparent performance.

**Feature set combinations.** Models are cross-validated over all 31 non-empty
subsets of the five features. Two further features (`readent`, `mhl`) were
dropped before the results were read, because their fill rate at the primary
depth was under 50 % — half of each column would have been other samples' means
rather than a measurement. That decision was made on fill rate, not performance,
and it reduced the combination count from 127 to 31.

---

## Input data

Not included, for size and access-control reasons. Point the two roots in
`0_setup/config.conf` at your copies.

**The cohort definition is also not included.** It lists sample identifiers, and
for controlled-access datasets the sample list may fall under the data access
agreement. A structure-only example is provided:

```
3_validate/steps/cohorts/sample_cohort.conf.example
  -> copy to ~/cohorts/<cohort>.conf and fill in the values
  -> add the cohort name to VAL_COHORTS in 0_setup/config.conf
```

| What | Where | Note |
|---|---|---|
| Training normals (BAM) | `<data root>/for_in_silico_test/normal_cfDNA_public/aligned_bam` | public cfDNA |
| Training cancer (cov, BAM) | `<data root>/training/GBM_cell-line/` | GBM cell lines |
| Validation cancer BAM | `CANCER` in the cohort definition | EGAD00001003427 (controlled access) |
| Validation normal background | `NORMAL` / `NORMALS_*` in the same file | zero overlap with training |

**The validation cancer samples are tissue DNA, not plasma cfDNA.** This is a
*simulation of clinical validation*, not clinical validation — confounding by the
tissue/cfDNA difference cannot be excluded. See the paper's limitations.

### What a fresh clone still needs

These inputs are required but **no script in this repository creates them**. The
error you get without each one is in the last column.

| Needed | Where the code looks | If missing |
|---|---|---|
| Normal cov files | `$METH_ROOT/sample_data/_input/normal_pub15/*.cov.gz` — under `METH_ROOT` (the *output* root), not `DATASET_ROOT`. Override with `NORMAL_SET=` | step 01: `정상 cov 가 없다` (both folders empty: `cov 파일 없음`) |
| Bisulfite genome index | `$METH_ROOT/data/Bisulfite_Genome/` | step 03 runs bismark **without** `--genome_folder`, which changes its behaviour silently |

Step `10c` builds the three panel-union CSVs and the two sample lists that step
11 consumes, from the panel files and `samples.csv`. Until 2026-10-01 those five
files were made by hand, so a fresh run stopped at step 11 with `패널 목록이
없다`. On the GBM data four of the five come out byte-identical to the hand-made
ones; `panel_union_rand3.csv` holds the same 586 rows in a different order, which
does not reach the output because `11_pool.py` ends in `samtools sort`.

BAM and cov filenames must share a sample id up to the first `_`. Step 03 pairs
them that way and refuses to run if any pair is missing (see below).

**Moving the data.** Changing `METH_ROOT` and `DATASET_ROOT` in `0_setup/config.conf`
moves everything: all 19 path blocks derive from those two. If your layout *under*
those roots differs, four environment variables move the data folders without
editing any file. Each value is written in two or three steps, so setting the
variable is safer than editing them by hand — edit two of three and the third
silently reads the old location.

| Variable | Default | Read by |
|---|---|---|
| `NORMAL_SET` | `$METH_ROOT/sample_data/_input/normal_pub15` | steps 01, 02, 03 |
| `GBM_COV_DIR` | `$DATASET_ROOT/training/GBM_cell-line/cov` | steps 01, 02, 03 |
| `NORMAL_BAM_DIR` | `$DATASET_ROOT/for_in_silico_test/normal_cfDNA_public/aligned_bam` | steps 03, 11 |
| `GBM_BAM_DIR` | `$DATASET_ROOT/training/GBM_cell-line/aligned_bam` | steps 03, 11 |

Step 03 pairs a cov file with its BAM by the sample id up to the first `_`. It
checks every pair **before** it starts scanning and stops if any is missing,
naming them and both BAM folders:

```
  GBM      짝 30개 · 못 찾음 0개
  Normal   짝 0개 · 못 찾음 15개
중단: 짝 BAM 을 못 찾은 검체가 15개 있습니다.
```

Any missing pair is fatal, not just an empty class, because step 01 writes
`samples.csv` from the cov list while step 03 scans only what it could pair.
Step 04 derives its normal-coverage threshold from `samples.csv`, so a cov set
of 28 against 15 BAMs either kills step 04 with `쓸 수 있는 창이 0` or, for a
smaller gap, quietly selects a different 200-block panel against a shrunken
normal reference. `ALLOW_SKIP=1` proceeds anyway and says so.

`--sample <id>` runs one sample and is exempt from the class accounting.

**Naming rules for a different dataset.** Two environment variables adapt the
sample-id and replicate-grouping rules so you never have to rename data files.
Both default to the conventions of the dataset used here, so leaving them unset
reproduces the reported results.

| Variable | Default | What it does |
|---|---|---|
| `SID_RE` | first `_`-separated field | Sample id from the filename. The first capture group is the id; a name the pattern does not match falls back to the default. |
| `GROUP_RE` | `SNU<n>-<rep>` → `SNU<n>` | Groups replicate files into one independent unit, **cancer side only**. Coverage requires a fraction of *units*, not of files, so this changes the threshold. |
| `PERSON_RE` | `pat_1_…` / `ind_3_…` | Stage 4 only. Groups validation samples by person so one person cannot carry several rows of a table. Prints a warning if it matches nothing. |
| `DEPTHS` | `5k 15k 50k 100k` | Panel depths to read in stages 3 and 4. Stage 3 also falls back to `_<n>k_` in the model name. |
| `CANCER_LABEL` | `GBM` | The value stage 1 writes in the `type` column of `samples.csv`, and the value every later step compares against. Changing it means re-running stage 1; a mismatch stops the run rather than training on zero positives. |

Steps 01 and 02 stop if two files resolve to the same sample id, because the id
becomes a column name and a collision would drop one of them silently.
`GROUP_RE` is not applied to normals: one file per person is assumed there and
an assertion enforces it.

For example, a dataset that separates sample from run with a double underscore
(`H1876_0d__SRR25758566_…`, `H1876__SRR25758603_…` — the same cell line
sequenced twice) runs with:

```bash
SID_RE='(.+?)__' GROUP_RE='(.+?)_\d+d$' bash 1_dmr/run_dmr.sh
```

giving distinct ids `H1876_0d` and `H1876` that group into one cell line.
A name the pattern does not match stays its own unit, so `SHP_77` is not
truncated to `SHP`. `GROUP_RE='([A-Za-z0-9]+)'` would truncate it.

Stages 2 and 4 take the normal/cancer split from the `type` column of
`samples.csv`, which step 01 writes. Nothing downstream infers class from a
filename. The cancer label itself is cosmetic — set `CANCER_LABEL=SCLC` and the
column reads `SCLC` — but it is not cosmetic in `2_train`, where sample ids are
built as `Normal1…` / `GBM1…` and the classifier keys on that prefix. Those ids
are internal and stay as they are.

**There are three pooled-BAM sets, not one.** A pooled BAM holds only reads inside
its own panel regions, so a null panel — whose blocks are drawn at random — finds
almost nothing in the real-panel pool. Stage 1 step 11 builds all three, and
`2_train/steps/config.py` picks the matching one from `PANEL` (`randb*` →
`_randb`, `rand*` → `_rand`, otherwise `_cellline`).

**Memory.** Stage 3 refuses to start the mixing step until `free -g` reports at
least 22 GB *available*, and gives up after about 3 hours. Only that step is
heavy: it uses roughly 1.35 GB per sample while the rest of the pipeline uses
almost none, which is why the gate sits there rather than on the driver. At most
`MEMGATE_N` samples (default 22) hold a token at once. Lower both with
`MEMGATE_FREE=<GB>` and `MEMGATE_N=<n>`; they are environment variables, not
config keys.

---

## Environment

Versions actually used are pinned in `requirements.txt`; `0_setup/ENVIRONMENT.md`
explains why each matters.

**Linux only.** The pipeline uses `flock`, `/proc`, `bismark`, `bowtie2` and
`samtools`. On Windows you can read the code and edit paths, but not run it.

```bash
conda create -n meth python=3.9
conda activate meth
pip install -r requirements.txt
# bismark, bowtie2, samtools are installed separately: see ENVIRONMENT.md
```

---

## Reproducibility

Re-running this code against the original inputs reproduced the stored outputs
byte for byte:

| Stage | What was re-run | Result |
|---|---|---|
| 1 | all nine panels, from candidate blocks | `DMR_confirmed_*.csv` md5 **9 / 9** |
| 2 | `bl200`, `jsd200`, `jsdb200` at 5k | `BEST.joblib` md5 **31 / 31** per panel |
| 3 | `bl200` at 5k, 178 samples | feature CSVs **6,230 / 6,230**, `y_prob_all_*.csv` **178 / 178** |
| 3 | `jsd200`, `jsdb200` at 5k, 4 samples each | feature CSVs **140 / 140**, `y_prob_all_*.csv` **4 / 4** |

Two caveats on the stage-2 comparison. The cross-validation tables
(`cv_*.csv`) for `jsd200` and `jsdb200` in our stored tree predate the decision to
drop `readent` and `mhl`, so they list 127 combinations where the re-run lists 31;
the models built from them are identical. And `04_coverage_check.py` writes a
diagnostic CSV whose row order used to follow directory-listing order; it now
sorts first, so that one file differs from copies written before the change.
Nothing reads it.

**A fixed seed is not a reproducibility guarantee.** Numbers can move if any of
these change, which is why `requirements.txt` pins what was used:

- library versions, especially `scikit-learn` (RandomForest and SVC defaults shift)
- BLAS thread count and `n_jobs` (floating-point accumulation order)
- `bismark` / `bowtie2` / `samtools` versions
- Python and OS

**Not re-run at full size.** `01_candidates.py` and `07_panel.py` were checked by
reading the code rather than by regenerating their inputs: `01_candidates` draws
no random numbers, and the models use `random_state=<bootstrap iteration>`.

**Small samples give different values.** `JSD_LIMIT` caps how many reads are
scanned. It is useful for checking that a stage runs, but the values will not
match; do not use it for value comparison.

---

## Things to know when reading the code

**Config lives in one place.** Edit `0_setup/config.conf` and every stage follows.
Environment variables win over the file — use them for isolated runs. Each output
folder gets a `config_snapshot.txt` recording, per key, where the value came from
(env / conf / default).

**File names are coupled to the config.** Validation *copies* the training scripts
into each sample's working folder. If you rename anything under `2_train/steps/`,
update `STEP_MAP` in `0_setup/config.conf` too. If you don't, that feature is
silently missing, scoring fills the gap with the training mean, and you get
**AUC 0.500**. A guard in the copy loop stops this.

**Generation (`GEN`)** is the suffix on panel output folders (currently `15`).
Mixing generations would let training and validation see different panels, so
several guards check it and stop on mismatch.

**Scripts refuse to be imported.** Calling an unguarded script through `importlib`
executes its body and overwrites real data. Run them as `python <file>` only.

**Stage 1 refuses to write into a non-empty output folder.** `1_dmr/run_dmr.sh`
checks, before running anything, every folder the selected steps create. Steps
that *add* to a folder an earlier step made (02, 04b, 04c, 09, 10) are exempt, so
they can be re-run alone. Stages 2–4 have no such guard: re-running
`2_train/run_train.sh` with the same `VERSION` overwrites the previous training
run without asking. To overwrite in stage 1 you must retype the folder name
(`OVERWRITE=<folder>`), which cannot be done by reflex.

**Locks and logs** live under `TMP_ROOT` (default `~/tmp`). Each driver creates
its own lock directory.

**Output names are ASCII.** Folders, files and identifiers are English throughout,
including everything the pipeline writes. `0_setup/RENAMES.txt` maps the script
names used in the paper's methods section to the names here.

> Code comments are in Korean and explain why each choice was made. The English
> README above covers what the code does.

---

## Layout

47 files. The four you run are marked ★.

```
README.md  requirements.txt  .gitignore  .gitattributes
0_setup/                    config.conf  <- the only file you edit
                            meth_config.py  IO.md  ENVIRONMENT.md  RENAMES.txt
1_dmr/     ★ run_dmr.sh     steps/  01_candidates  02_blockcpg  03_jsdcount  04_jsd
                                    04b_jsdnull  04c_jsdboot  05_llrref
                                    06_select  06b_select_supplement  07_panel
                                    08_export_panels  09_null_a  10_null_b
                                    11_pool  run_pool.sh  moderated_t
2_train/   ★ run_train.sh   steps/  01-03_Preprocessing…Mixing  04_coverage_check
                                    05_Feature_Matrix  06_jsdfeat  07_llrfeat
                                    08_ML_classifier  09_save_best_j  config
3_validate/★ run_validate.sh steps/ run_one_sample  make_val_config  score_one_sample
                                    paths  val_single  cohorts/
4_report/  ★ make_report.sh steps/  primary_metric  make_tables
```
