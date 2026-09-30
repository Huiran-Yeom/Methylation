# Methylation-based GBM classifier: panel selection to external validation

Code for building a methylation-pattern classifier for glioblastoma (GBM) from
cfDNA and validating it on an external cohort. This folder contains only the code
that produces the numbers reported in the paper.

> **Quick check.** `bash selfcheck.sh` prints a PASS/FAIL table over 12 groups and
> writes nothing. It currently reports **PASS 26 · FAIL 0**. The groups cover
> structure, step counts, syntax of every `.py` and `.sh`, config resolution, naming,
> cross-references, publishable identifiers, heredoc exports, output paths,
> the step-name map, stray files, and where paths are declared.
>
> Groups 6, 8, 9 and 10 exist because each caught a defect that had already shipped:
> a driver calling a renamed file, a variable unexported into a here-document, a
> Korean path fragment, and a step-name map that had silently stopped being read.
> Each was re-broken on purpose to confirm the group turns red.

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
| **1** | DMR panel selection | `1_dmr/run_dmr.sh` | candidate blocks → JSD → 200-block panel + 6 null panels |
| **2** | Training | `2_train/run_train.sh` | in-silico mixing → five features → four models, cross-validated |
| **3** | External validation | `3_validate/run_validate.sh` | mix the external cohort under identical conditions, then score |
| **4** | Reporting | `4_report/make_report.sh` | discrimination tables · paired bootstrap intervals |

Each stage folder holds exactly **one entry point**; its step code sits in `steps/`,
numbered in dependency order. There are more step files than stages because these
are the files that produced the reported numbers — they were not merged or rewritten
for publication, so the md5 checks in *Reproducibility* below refer to this exact
code. `0_setup/IO.md` lists every step with its inputs and outputs if you need to
go deeper; you do not need it to run the pipeline.

**Stage 4 reports nine panels at every depth; stages 2 and 3 run one panel at one
depth.** Stage 1 builds all nine panels (3 real + 6 null) in a single run and does
not depend on depth. Training and validation work on whichever `VERSION`/`PANEL`
you point them at, so the full report needs stages 2–3 run once per panel *and*
per depth. The reported tables are **15 combinations**: the three real panels at
5k / 50k / 100k, and the six null panels at 5k only.

Depth lives in two places that must agree — the `5k` inside `VERSION` and
`TARGET_DEPTH_READS`. They are checked: `0_setup/meth_config.py` stops if they
disagree, and `make_val_config.py` stops if `VERSION` and the panel you pass to
stage 3 disagree. Environment variables win over `config.conf`, so the loop below
does not edit any file:

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

That is 15 training runs and 15 validation runs. On our server one validation
panel-depth took about 33 minutes per sample, so 178 samples at 40-way parallelism
is roughly 2 hours per combination.

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
| Null B `randb1-3` | stratified draw from the Baseline candidate pool, **matching the Baseline CpG-count distribution** | Baseline |

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
`MUT_RATIOS` are 0, 0.1, 0.5, 1, 2, 2.5, 3, 5, 10 and 100%. The RNG
is seeded from `SEED`, the ratio, the replicate **and the sample name**, so each
combination reproduces and no two samples share a negative draw. (The sample axis
was added 2026-09-16 after all 178 validation negatives came out read-identical.) Ratios whose pool is too shallow are **skipped rather than filled by
sampling with replacement**: filling would make the replicates non-independent
and inflate apparent performance.

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

### What a fresh clone still needs

An adversarial review walked the five commands on a clean machine. These are the
inputs the pipeline requires that **no script in this repository creates**. They
existed on our machine from earlier work; if you start from nothing you must
supply them, and the error you get otherwise is named in the last column.

| Needed | Where the code looks | If missing |
|---|---|---|
| Normal cov files | `$METH_ROOT/sample_data/_input/normal_pub15/*.cov.gz` — note this is under `METH_ROOT` (the *output* root), not `DATASET_ROOT`. Override with `NORMAL_SET=` | step 01: `cov 파일 없음` |
| **Three** panel-union CSVs | `$METH_ROOT/results/dmr/j_panel<GEN>/` — one `chr,blk` row per panel block, assembled by hand from stage-1 output. `panel_union3_cellline.csv` (the three real panels), `panel_union_rand3.csv` (null A), `panel_union_randb3.csv` (null B). Each feeds one `run_pool.sh` mode. | step 11: `패널 목록이 없다` |
| `use<GEN>_normal.txt`, `use<GEN>_gbm.txt` | `~/tmp/` — one sample id per line (the id is the filename up to the first `_`). Override with `USE_NORMAL=` / `USE_GBM=` | step 11: `쓸 BAM 목록이 없다` |
| Bisulfite genome index | `$METH_ROOT/data/Bisulfite_Genome/` | step 03 runs bismark **without** `--genome_folder`, which changes its behaviour silently |

BAM and cov filenames must share a sample id up to the first `_`; step 03 matches
them that way and **skips** any cov whose BAM it cannot find.

**Memory.** `3_validate` refuses to start the mixing step until `free -g` reports
at least 22 GB *available*, and gives up after ~3 hours. On a smaller machine every
sample fails after that wait. Lower it with `MEMGATE_FREE=<GB>` and `MEMGATE_N=<n>`
(concurrent slots): both are environment variables, not config keys.

| What | Where | Note |
|---|---|---|
| Training normals (BAM) | `<data root>/for_in_silico_test/normal_cfDNA_public/aligned_bam` | public cfDNA |
| Training cancer (cov, BAM) | `<data root>/training/GBM_cell-line/` | GBM cell lines |
| Validation cancer BAM | `CANCER` in the cohort definition | EGAD00001003427 (controlled access) |
| Validation normal background | `NORMAL` / `NORMALS_*` in the same file | zero overlap with training |

**The validation cancer samples are tissue DNA, not plasma cfDNA.** This is a
*simulation of clinical validation*, not clinical validation — confounding by the
tissue/cfDNA difference cannot be excluded. See the paper's limitations.

---

## Environment

Versions actually used are pinned in `requirements.txt`; `0_setup/ENVIRONMENT.md`
explains why each matters.

**Linux only.** The pipeline uses `flock`, `/proc`, `bismark`, `bowtie2`, and
`samtools`. On Windows you can read the code and edit paths, but not run it.

```bash
conda create -n meth python=3.9
conda activate meth
pip install -r requirements.txt
# bismark, bowtie2, samtools are installed separately — see ENVIRONMENT.md
```

---

## Reproducibility: what is and is not guaranteed

Re-running with the same input and config produced byte-identical output
(one sample, one panel, verified repeatedly):

| Target | Result |
|---|---|
| 9 panels (`DMR_confirmed_*.csv`) | md5 9/9 |
| Feature CSVs, one sample × 8 ratios (mean, entropy, jsd, pdr, llr) | md5 **39/40** |
| — the 40th (`coverage_check_*.csv`) | content identical, row order differs |
| Per-sample scores (`y_prob_all_*.csv`) | md5 identical (4,340 rows) |
| Tables (9 panels × 2 cohorts; 3 depths) | output identical |
| Config provenance in `config_snapshot.txt` | every key resolved from `conf` |

The one row-order difference is deliberate. `04_coverage_check.py` used to sort
ties by insertion order, which follows directory-listing order, so its diagnostic
CSV was not byte-reproducible. It now sorts the index first and stable-sorts by
count, which fixes the order but makes it differ from files written before the fix.
Nothing reads that CSV; the five feature matrices are unaffected.

**A fixed seed is not a reproducibility guarantee.** Numbers can move if any of
the following change — which is why `requirements.txt` pins what was used:

- library versions (especially `scikit-learn` — RandomForest and SVC defaults shift)
- BLAS thread count and `n_jobs` (floating-point accumulation order)
- `bismark` / `bowtie2` / `samtools` versions
- Python and OS

**Not directly verified**

- `1_dmr/steps/01_candidates.py` and `07_panel.py` were **not re-run at full input
  size.** Seeds were confirmed by reading the code (`RandomState(42)`; models use
  `random_state=<bootstrap iteration>`; `01_candidates` draws no random numbers).
  So "the inputs did not change" is verified; "regenerating the inputs gives the
  same thing" is not.
- `11_pool.py` (panel-region pooled BAMs) was **not executed by us** — the pooled
  BAMs it writes were produced once, by hand, before this driver existed. Stage 2
  step 01 reads them (`config.py` → `FULL_GBM_BAM` / `FULL_NORMAL_BAM`), so it is
  on the critical path, not optional. Its inputs, the three `panel_union*.csv`
  files, are **not written by any script here**: ours were assembled by hand from
  the stage-1 panel outputs. If you start from scratch you must create all three
  (one `chr,blk` row per panel block) before stage 11 will run.
- **There are three pooled-BAM sets, not one.** A pooled BAM holds only reads
  inside its panel regions, so a null panel — whose blocks are drawn at random —
  finds almost nothing in the real-panel pool. Measured 2026-09-30 when we tried
  it: `rand1200` passed **8 of 200** blocks, `rand3200` stopped with
  `풀 부족으로 중단`. `config.py` therefore picks the pool from `PANEL`
  (`randb*` → `_randb`, `rand*` → `_rand`, otherwise `_cellline`), matching the
  per-panel configs the reported runs used. Build all three:
  `bash 1_dmr/steps/run_pool.sh cellline` / `rand` / `randb`.
- `04b_jsdnull.py` / `04c_jsdboot.py` / `06b_select_supplement.py` produce the
  `jsdb200` panel. They ran once, on 2026-09-10, before the numbered layout existed;
  the panel CSV they wrote is what the reported `jsdb200` results come from, and
  it was **not regenerated** during verification.

**Small samples give different values.** `JSD_LIMIT` caps how many reads are
scanned, which is useful for checking that a stage runs — but the values will not
match. Do not use it for value comparison.

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
**AUC 0.500**. A guard in the copy loop now stops this.

**Generation (`GEN`)** is the suffix on panel output folders (currently `15`).
Mixing generations would let training and validation see different panels, so
several guards check it and stop on mismatch.

**Scripts refuse to be imported.** Calling an unguarded script via `importlib`
executes its body and overwrites real data. Run them as `python <file>` only.

**Locks and logs live under `~/tmp/lock/` and `~/tmp/`.** Each driver `mkdir -p`s
its own lock directory. In our internal tree that directory happened to be created
by a monitoring script which is not published here, so a fresh clone used to fail
on the first `flock`; the drivers now create it themselves.

**Stage 1 refuses to write into a non-empty output folder.** `1_dmr/run_dmr.sh`
checks, before running anything, every folder the selected steps write to.
Stages 2–4 have no such guard: re-running `2_train/run_train.sh` with the same
`VERSION` overwrites the previous training run without asking. To overwrite you must
retype the folder name (`OVERWRITE=<folder>`), which cannot be done by reflex.

**Output names are all ASCII.** Folders, files, and identifiers are English
throughout, including the names of everything the pipeline writes. Our internal
tree used Korean output names while the analysis was run; writers and readers were
changed together, and `selfcheck.sh` group 9 checks that none are left.
`0_setup/RENAMES.txt` gives the old↔new table so the paper's methods can be traced.

> Code comments are in Korean. They record *why* each choice was made — including
> several cases where an earlier version was silently wrong — and translating them
> would lose that. The English README above covers what the code does.

---

## Layout

47 files. The four you run are marked ★.

```
selfcheck.sh                one command: PASS/FAIL table (read-only)
README.md  requirements.txt  .gitignore
0_setup/                    config.conf  <- the only file you edit
                            meth_config.py  IO.md  ENVIRONMENT.md  RENAMES.txt
1_dmr/     ★ run_dmr.sh     steps/  01_candidates  02_blockcpg  03_jsdcount  04_jsd
                                    04b_jsdnull  04c_jsdboot  05_llrref
                                    06_select  06b_select_supplement  07_panel
                                    08_export_panels  09_null_a  10_null_b
                                    11_pool  run_pool  moderated_t
2_train/   ★ run_train.sh   steps/  01-03_Preprocessing…Mixing  04_coverage_check
                                    05_Feature_Matrix  06_jsdfeat  07_llrfeat
                                    08_ML_classifier  09_save_best_j  config.py
3_validate/★ run_validate.sh steps/ run_one_sample  make_val_config  score_one_sample
                                    paths  val_single  cohorts/<example>
4_report/  ★ make_report.sh steps/  primary_metric  make_tables
```

`RENAMES.txt` maps every file back to the name it had while the analysis was run, so
the paper's methods section can be traced to this tree, and lists what is **not**
here: operations and monitoring scripts, and the parallel/resume launchers. Those
do not touch a reported number. Everything that does is in this tree, including
stage 11 and the two JSD null/bootstrap steps.

Operations and monitoring scripts (memory watchdog, progress board, resume
drivers) are not included: they do not affect results.
