# Unified CRISPR Off-Target Workflow

A self-contained arm added to nf-core-scge on branch `feat/offtarget-wgs`. It runs
via a **named entry** (`-entry OFFTARGET`) and **does not touch the default SCGE pipeline**.

(2026-07-09): the legacy `crispr_ml_*` read-level classifier is **deprecated** for
off-target work; this workflow uses the pileup shape-model approach
(`worklist_from_vcf` / `pon_filter` / `score` + `wgs_shape_model.pkl`).

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="images/offtarget_metro_dark.svg">
  <img alt="Unified CRISPR Off-Target Workflow metro map" src="images/offtarget_metro.svg" width="820">
</picture>

*Three lines: **ECS truth** (red) at hotspots, **WGS genome-wide discovery** (blue) → PoN-filtered
worklist, and the **WGS-only model** (green) that joins WGS features to ECS truth. Source:
[`offtarget_metro.mmd`](offtarget_metro.mmd) (rendered with [nf-metro](https://github.com/seqeralabs/nf-metro);
an interactive pan/zoom version is at [`images/offtarget_metro.html`](images/offtarget_metro.html)).*

## What you get

Everything lands in `<outdir>/offtarget/`.

- **`wgs_offtarget_worklist_pon.csv`** — the main result. Every candidate edit found in the WGS,
  ranked, with germline and artifact calls already knocked out by the panel of normals. Start here.
- **`offtarget_report.csv`** — the same worklist with two extra columns: whether each hit falls on a
  known hotspot, and whether ECS confirmed it.
- **`<id>.offtarget_analysis.tsv`** — the ECS answer at each hotspot: edited or not, and at what VAF.
  This is the ground truth.
- **`recall_vs_vaf.csv` / `.png`** — how often the WGS catches an ECS-confirmed edit, split by VAF.
  Read it as: above this VAF, believe the WGS calls; below it, don't.
- **`offtarget_metrics.json` / `.txt`** — PR-AUC and recall-weighted F-beta (F2, F5) for the WGS
  shape score, plus precision/recall at the reported operating point. `.txt` is the one-screen
  human read; `.json` is the machine-readable copy. **Read the denominator section first** —
  the precision in here is against the ECS label, *not* against human review.
- **`training.tsv`** — one row per hotspot lining up the WGS features against the ECS answer, with a
  `label` marking genuine *somatic* edits: ECS shows an edit **and** it's absent from the matched WGS
  normal (germline/recurrent sites are demoted, so they don't get called edits). Only used to retrain
  the model offline; ignore it on a normal run.

The last three only show up when you give it both ECS and WGS.

### Verification snapshots (optional)

Add `--offtarget_snapshots true` and every **LIKELY EDIT** in the WGS worklist gets an IGV-style
read pileup — **edited (tumor) on the left, matched normal on the right** — written to
`<outdir>/offtarget/snapshots/`, one PNG per call. A real edit shows indel-bearing reads (red
deletions, purple insertion ticks) stacked at the locus in the tumor over a clean normal; an
artifact shows up in both or neither. It's the by-eye confirmation step, straight from the CRAM,
no IGV needed. Off by default (`offtarget_snapshots = false`) because it renders one image per hit.

![tumor vs normal pileup snapshot](images/offtarget_snapshot_example.png)

*Example output: the AAVS1 on-target (chr19:55,115,731). Left (edited) — deletions and insertions
pile up at the cut site; right (unedited normal) — clean. Pass `--snapshots` to `run_offtarget.sh`
to turn this on, so the confirmed off-targets (PLCB2, CNNM3) get the same tumor/normal packet.*

### Read-level tags for IGV (optional)

`<id>.offtarget_analysis.tsv` gives you `indel_fraction = 0.42` and leaves you to work out by eye
which reads made up the 0.42. Add `--offtarget_tagged_bam true` and the ECS caller also writes
`<id>.tagged.bam` (+ `.bai`), in which **every read carries an `XC` string tag naming how the caller
classified it** — so the pileup shows you its reasoning directly.

In IGV: load the BAM → right-click the track → **Color alignments by** → **tag** → type `XC`.

| tag | meaning |
|---|---|
| `Edited_Deletion_<N>bp` | deletion of N bp called from the CIGAR or a split alignment |
| `Edited_Insertion_<N>bp` | insertion of N bp |
| `Edited_Duplication_<N>bp` | tandem duplication from a split alignment |
| `Edited_BND_<chrom>` | breakend — the read's other end maps to `<chrom>` (e.g. the transgene contig) |
| `Edited_SoftClip` | event recovered by realigning a soft clip |
| `Unedited_WT` | spans the target cleanly, no event — this is the denominator |
| `Skipped_Duplicate` / `_LowMapQ` / `_Mismatches` / `_Secondary` / `_Supplementary` / `_Unmapped` | excluded by a read filter, shown so you can see *why* it isn't counted |
| `Skipped_Unevaluable` | a classifier ran but found nothing near enough to a PAM to call |
| `Skipped_NoSpan` | inside the padded window but doesn't span the target and carries no event |

Each alignment appears **exactly once**, so the depth IGV shows is real. Where a read falls in two
overlapping target windows and gets classified differently, the most specific call wins — the table
above is in precedence order, top to bottom.

**Don't expect the tag counts to equal the TSV columns** — they count different things, and the BAM
is the more literal of the two:

- tags are per **alignment record**, while `indel_reads` is per **fragment** (the caller collapses
  R1/R2 by read name) and is taken *after* the site-level filters. So `Edited_*` records run higher
  than `indel_reads` — roughly 2× where both mates cover the cut site. On the AAVS1 site14
  on-target that is 17,396 `Edited_*` records against `indel_reads` = 10,821.
- the BAM covers the whole padded window, so it also holds reads that never entered `total_reads`.
  Across the AAVS1 site14 panel that is 65% `Unedited_WT`, 35% `Skipped_*` and 0.7% `Edited_*`; at
  a *1 bp* target specifically, the ±150 bp pad means most records are `Skipped_NoSpan` (69% at the
  on-target). They are kept deliberately — a pileup cropped to reads spanning a single base is
  unreadable in IGV.

Use the tags to see *which* reads drove a call and why; use the TSV for the number.

**Off by default, and worth keeping that way for routine runs.** Output is restricted to the target
windows (target ±150 bp) rather than the whole CRAM, but it still scales with target count and
depth. Measured on one AAVS1 site14 ECS sample (1,149 target intervals → 1,145 merged windows,
~11,000× at the on-target):

| | without | with `--offtarget_tagged_bam` |
|---|---|---|
| output | — | **800 MB** `.bam` + 1.5 MB `.bai`, 21.8 M reads |
| peak RSS | 1.4 GB | **6.3 GB** (the tag map is held until the write pass) |
| wall clock | ~30 min | ~46 min |

`ECS_INDELS` is given 24 GB instead of 8 GB when the flag is set, so enabling it does not OOM the
truth arm. Still check the size on one sample before turning it on across a cohort — and note the
debug log this pipeline learned that lesson from, `offtarget_ecs_unevaluable_log`, reached
0.1–1 TB/sample and filled the work directory. Restricting the run with `--regions` keeps both
numbers small when you only want to review a handful of loci.

## What runs depends on the samplesheet

The `datatype` column decides:

- **ECS + WGS rows** → the whole thing: worklist, report, training table, recall curve. Use this to
  build or check the model.
- **WGS only** → worklist and report. No training table or recall curve.
- **ECS only** → just the hotspot truth tables.

The run logs the mode it picked up front (`OFFTARGET mode: paired / wgs_only / ecs_only`), and it
stops early with a clear error if the `datatype` column is missing or has anything other than
`ecs`/`wgs`. For **WGS rows** it also preflights the DRAGEN sidecars: the matched normal
(`<base>.cram`) and somatic VCF (`<base>.hard-filtered.vcf.gz`) are derived from the tumor CRAM
name (`<base>_tumor.cram`) and must sit beside it — if any is misnamed or missing, the run fails
immediately with the exact file and row, instead of a confusing empty worklist or a mid-run crash.

## Running it

On RIS Compute2 (SLURM + Apptainer) — the validated path. One wrapper handles every cohort; pass
the samplesheet and an output dir (`--help` lists the options):

```bash
sbatch run_offtarget.sh --input <samplesheet.csv> --outdir <dir> [--snapshots]
# the head job runs here and submits the task jobs to SLURM
```

or directly (the head job must run somewhere that can `sbatch` — a login/compute node, **not** the
interactive JupyterLab exec node):

```bash
module load nextflow apptainer
nextflow run . -entry OFFTARGET -profile ris2,apptainer \
    --input offtarget_samplesheet.csv \
    --outdir ./results_offtarget -resume
```

On RIS Compute1 (LSF) run `nextflow run . -entry OFFTARGET -profile ris` under `bsub` instead.
There's a filled-in example samplesheet at `assets/offtarget_samplesheet_template.csv`.

You hand it one samplesheet with these columns:
`sample,datatype,guide,edited_cram,control_cram,target_file,vcf`. A few rules:

- `datatype` is `ecs` or `wgs`. That's what sorts each row into the right arm.
- `guide` links an ECS sample to its WGS partner — same guide means same experiment.
- For a WGS row, point `edited_cram` at the DRAGEN `<name>_tumor.cram`, and keep `<name>.cram` (the
  normal) and `<name>.hard-filtered.vcf.gz` in the same folder. The scripts find them by name.
- Use absolute paths — the CRAMs are read straight off storage, not copied in.

## Key params (in `nextflow.config`)

| param | default | purpose |
|---|---|---|
| `offtarget_shape_model` | `assets/models/wgs_shape_model.pkl` | pileup shape ranker (fixed asset; train offline) |
| `offtarget_min_af` | 0.05 | DRAGEN tumor AF floor for candidate indels |
| `offtarget_min_span` | 8 | spanning-read coverage gate |
| `offtarget_hi_score` | 0.60 | shape score ≥ this = "detected" for the recall curve |
| `offtarget_target_recall` | 0.80 | recall level whose VAF floor is reported |
| `offtarget_germline_max_ctrl_if` | 0.05 | matched-normal indel frac above this = germline/artifact, not a somatic edit (`label` 0) |
| `offtarget_hotspot_pad` | 25 | bp window to match a worklist hit to a predicted hotspot |
| `offtarget_snapshots` | false | render IGV-style pileup PNGs for LIKELY EDITs |
| `offtarget_tagged_bam` | false | emit `<id>.tagged.bam` with per-read `XC` tags for IGV review (above) |
| `offtarget_rescue` | true | high-evidence rescue (below); `false` = call on model score alone |
| `offtarget_rescue_min_ifrac` | 0.15 | rescue: minimum indel fraction in the edited sample |
| `offtarget_rescue_min_conc` | 0.5 | rescue: minimum positional concordance (clonality) |
| `offtarget_rescue_min_span` | 20 | rescue: minimum spanning reads |
| `offtarget_metrics_betas` | `2,5` | F-beta weights for `offtarget_metrics.{json,txt}` |
| `offtarget_metrics_negatives` | `ecs_negative` | PR-AUC negative set; `all_label0` also counts germline-demoted ECS edits |

## Reading `recall_vs_vaf.csv` (the denominator matters more than the number)

A recall figure is only as honest as the set of "real edits" it divides by. The raw ECS label
counts **any** nonzero indel fraction as an edit (`offtarget_ecs_edit_threshold` defaults to
`0.0`), and at ECS depth that is overwhelmingly noise — in a real AAVS1 run, 10,780 of the
12,067 sites with any indel sat below 0.5% VAF with a **median of 3 indel reads out of ~5,000**.
Dividing by those produced a headline recall of ~0.002 that measured nothing except how much
ECS noise the pipeline correctly ignores.

The denominator is now a **credible ECS edit** — somatic (`label == 1`), `ecs_if >=
offtarget_min_ecs_vaf`, and `ecs_indel_reads >= offtarget_min_ecs_reads`. Read support is the
load-bearing half: VAF alone cannot tell a genuine 0.5% edit at 5,000× from 3 stray reads.

Sites with no spanning WGS reads are **unevaluable, not missed** — the scorer never got to see
them — so they are reported separately instead of being charged against recall:

| column | meaning |
|---|---|
| `n` | credible ECS edits in the bin |
| `n_unevaluable` | of those, sites with no spanning WGS reads (the depth floor) |
| `n_evaluable` / `n_detected` | sites WGS could judge / of those, called LIKELY EDIT |
| `recall` | `n_detected / n_evaluable` — scoring performance |
| `recall_incl_unevaluable` | `n_detected / n` — pessimistic view, depth floor included |
| `denom_min_ecs_vaf`, `denom_min_ecs_reads`, `denom_excluded_as_noise` | the denominator definition, stamped in so the file is self-describing |

On the AAVS1 run this turns a meaningless `0.002` into an interpretable curve: the top VAF bin
is 4 credible edits, 2 of them without WGS coverage, and **2/2 of the evaluable ones detected**.
Low-VAF bins still read low — that is the genuine WGS depth floor, and it is the honest result.

## Reading `offtarget_metrics.txt` (PR-AUC and F2/F5)

A recall-first diagnostic is judged on PR-AUC and recall-weighted F-beta, not on ROC-AUC. At the
positive prevalence this assay works at, ROC-AUC is dominated by the true-negative mass and reads
flatteringly high; PR-AUC has no such property. It is still reported, just not as the headline.

```
F_beta = (1 + beta^2) * P * R / (beta^2 * P + R)
```

beta weights recall `beta^2`× more than precision (β=2 → 4×, β=5 → 25×). F_beta is **monotone in
beta** — it rises with beta when recall > precision and falls when precision > recall — so F1 is
always an *endpoint* of {F1, F2, F5}, never the middle value. F1 in the middle means the beta
wiring is inverted.

### You cannot compute precision from the manual review

The human-reviewed table (`cart_ecs_merged.csv.gz`) has **55** rows with `manual_review == '1'`,
one `'1?'`, and **80,440 NaN**. `NaN` means *not reviewed* — not reviewed-and-rejected. There are
no confirmed negatives in it, so precision against it is not defined. Filling those NaNs with 0
would count every genuine discovery the reviewers never got to as a false positive and produce a
confidently wrong number. That is why `validate_recall.py` reports **recall only**, and it stays
that way.

So there are **two denominators**, kept strictly apart and both stamped into every output:

| | denominator | what it supports | where |
|---|---|---|---|
| 1 | ECS label in `training.tsv` (a real two-class label) | PR-AUC, precision, F-beta | `offtarget_metrics.{json,txt}` |
| 2 | human manual review (positives only) | **recall only** — the gold standard | `bin/validate_recall.py` |

The positive set uses the same credibility gates as the recall curve (`ecs_if >=
offtarget_min_ecs_vaf`, `ecs_indel_reads >= offtarget_min_ecs_reads`); without them the positives
are re-polluted with sub-0.5%-VAF ECS noise and every metric becomes meaningless. `label == 1`
rows that fail those gates are **ambiguous, not negative** — ECS saw something below the
credibility floor — so they are excluded and counted, for exactly the same reason the unreviewed
manual-review rows are. Rows with no score (`INSUFFICIENT COVERAGE` / `NO CRAM`) are excluded from
the ranking metrics and the exclusion count is reported; they are never filled with 0.

`label == 0` is likewise two populations: genuine ECS-negatives, and ECS edits **demoted because
the indel is in the matched normal**. The germline-demoted rows carry a real indel, so the shape
ranker correctly scores them high — germline rejection is a separate downstream gate, and none of
them are called `LIKELY EDIT`. Including them in the *ranking* denominator charges the ranker with
a job it does not do, so `offtarget_metrics_negatives` defaults to `ecs_negative`. The alternative
is always printed as a labelled sensitivity block, so the choice is visible rather than buried.

### What the AAVS1 cohort actually shows

```
n_pos / n_neg : 631 / 299     prevalence 0.679
PR-AUC        : 0.655  (0.96x prevalence)
operating pt  : TP 21  FP 5  FN 610  TN 294
                precision 0.808   recall 0.033   F1 0.064   F2 0.041   F5 0.035
```

PR-AUC **below** prevalence normally means a useless ranker or a bug. Here it is neither, and the
by-VAF-floor table in the same file shows why — 427 of the 631 credible positives sit in the
0.5–1% VAF bin, which WGS at ~30× genuinely cannot see:

| ECS VAF floor | n_pos | prevalence | PR-AUC | lift | ROC-AUC |
|---|---|---|---|---|---|
| ≥ 0.005 | 631 | 0.678 | 0.655 | 0.96× | 0.439 |
| ≥ 0.01 | 204 | 0.406 | 0.456 | 1.12× | 0.488 |
| ≥ 0.02 | 58 | 0.162 | 0.305 | 1.87× | 0.622 |
| ≥ 0.05 | 14 | 0.045 | 0.290 | **6.49×** | 0.931 |
| ≥ 0.10 | 9 | 0.029 | 0.257 | **8.81×** | 0.939 |

The ranker is strongly informative exactly where WGS can see (≥5% VAF: 6.5× lift, ROC 0.93) and
uninformative at the depth floor. That is the same conclusion the recall curve reaches — recall
`0.033` over 631 evaluable sites here reconciles exactly with `recall_vs_vaf.csv` — stated in
precision/recall terms. A **flat** lift at every VAF floor would be the real alarm.

Precision `0.808` is against the ECS label. It is **not** a claim about precision versus human
review, and the file says so on its own face.

## The high-evidence rescue (why recall is 100%, not 90%)

The shape model is a **ranker**, not a detector. Measured against the human-reviewed gold
standard (`Manual Indel Review/cart_ecs/cart_ecs_merged.csv.gz`, `manual_review == 1`), scoring
on the model alone recovered **47/52** confirmed CART edits. The five misses were not close
calls — indel fractions of 0.25–1.00, `ctrl_if` 0.00, and 130–244 spanning reads, i.e. pileups a
reviewer calls instantly in IGV — that the ranker scored 0.27–0.50. A sixth was blocked by the
hard `conc_ratio ≥ 0.5` clonality gate despite a 0.99 model score: several indel alleles sharing
one cut site, which is what multi-allelic Cas9 editing looks like.

So `score.py` calls **LIKELY EDIT** on unambiguous evidence regardless of model score, via two
arms, both gated *behind* the low-MAPQ and matched-normal vetoes — the rescue can override the
**model**, never the **evidence**, so it cannot resurrect a germline variant or a repeat pile-up:

1. **under-scored** — `indel_frac ≥ 0.15`, `conc_ratio ≥ 0.5`, `spanning ≥ 20`.
2. **multi-allelic** — model score ≥ `offtarget_hi_score` with the depth/burden evidence of (1),
   waiving the clonality requirement.

`call_basis` in `wgs_hotspot_scores.csv` records which fired (`model` / `high-evidence`).

**Measured cost.** Across both real cohorts the rescue adds calls only where they belong:

| cohort | scored rows | LIKELY EDIT before → after | recall vs manual review |
|---|---|---|---|
| CART | 99,308 | 91 → 97 (+5 under-scored, +1 multi-allelic) | 47/52 → **52/52 (1.000)** |
| AAVS1 | 5,088 | 34 → 35 (+1 under-scored) | on-target control preserved |

All six extra CART calls are the six confirmed edits; no other row in 99,308 was promoted.
Recall is 49/49 on-target and **3/3 off-target**.

Reproduce (and gate CI) with:

```bash
python bin/validate_recall.py \
    --scores results/offtarget/wgs_hotspot_scores.csv \
    --training results/offtarget/training.tsv \
    --gold "Manual Indel Review/cart_ecs/cart_ecs_merged.csv.gz" \
    --require-recall 1.0
```

Set `--no-rescue` (or `--offtarget_rescue false`) to reproduce pre-rescue behaviour.

> **Scope.** 100% is recall against the *human-reviewed* truth set — every edit a reviewer
> confirmed, the pipeline reports. It is not a claim that no edit exists below the WGS depth
> floor: an ECS edit with no spanning WGS reads is still unrecoverable and surfaces as
> INSUFFICIENT COVERAGE. The rescue thresholds were tuned on these 52 sites and re-checked on
> AAVS1; they should be re-validated against any new cohort's manual review.

## Retraining the shape model (`-entry TRAIN`)

The deployed model (`offtarget_shape_model`) is a fixed asset — training is deliberately **not**
in the OFFTARGET DAG. To fit a new one from your own cohort, use the separate `TRAIN` entry:

```bash
# 1. run OFFTARGET (paired ECS+WGS) to produce the labeled training table
sbatch run_offtarget.sh --input paired_samplesheet.csv --outdir results
#    -> results/offtarget/training.tsv

# 2. fit a new shape model from it
nextflow run . -entry TRAIN -profile ris2,apptainer \
    --input results/offtarget/training.tsv --outdir results
#    -> results/train/wgs_shape_model.pkl  +  train_metrics.json (holdout AUC/AP)

# 3. deploy it
sbatch run_offtarget.sh --input cohort.csv --outdir results2 \
    -- --offtarget_shape_model results/train/wgs_shape_model.pkl
```

`train_shape_model.py` fits a `HistGradientBoostingClassifier` on the `MODEL_FEATURES` **present in
the table** and records that exact list in the bundle (`{model, features, …}`); `score.py` selects
inputs by `bundle["features"]`, so the model is a drop-in even if the table carries a subset. It
runs in the off-target container so the pickle is written under the same scikit-learn the pipeline
scores with (the version guard rejects loading a model pickled under a *newer* sklearn). Tune with
`--offtarget_train_{learning_rate,max_iter,max_depth,holdout_frac,seed}`.

> The full-strength model needs all seven features (`indel_frac, conc_ratio, pos_conc, pos_mad,
> modal_len, modal_mapq, softclip_frac`); `SCORE_HOTSPOTS` now emits all of them into `training.tsv`.
> A table generated before that change carries only three and trains a weaker model (the trainer
> warns and reports which are missing).

## What this workflow is (and is not)

- **It is** a **hotspot edit-confirmation + genome-wide screen**. At known/nominated hotspots the WGS
  shape scorer recovers edits well; genome-wide it produces a PoN-filtered, ranked worklist for review.
  On-target recovery is proven end-to-end: on the first real AAVS1 run the AAVS1 on-target
  (chr19, *PPP1R12C*) came back with `is_hotspot=1, ecs_confirmed=1` for both guides, **from WGS alone**.
- **It also** recovers the one confirmed *off*-target edit we have — **PLCB2 chr12:32,679,410** (90% VAF,
  ECS-confirmed, ~202 WGS indel reads). That is the proof WGS-only can find a real, homology-based
  off-target. See `OFFTARGET_POSITIVE_CONTROL.md` at the project root.
- **It is not (yet)** a workflow with a *quantified* off-target detection floor ("trust down to X% VAF").
  Across both cohorts (AAVS1 + the 30-guide CAR-T screen) the human-reviewed truth set holds only
  **2 confirmed off-targets**, both >85% VAF — there is no low-abundance off-target population to measure
  recall against. That's the editing being highly on-target-specific, not a bug; but it means sub-5%
  off-target sensitivity is unproven. Treat WGS-only calls as trustworthy at hotspots **≥5% VAF**.

## Validation

**First real run** (AAVS1, 2 WGS × 6 ECS, RIS Compute2 / SLURM + Apptainer): completed end-to-end; the
ECS⋈WGS join produced a labeled `training.tsv`; the on-target was recovered in WGS alone (above). The
recall-vs-VAF curve behaves as depth predicts — near-zero below ~1% VAF, rising with abundance — but the
off-target bins above 5% are too sparse (single digits) to fix a floor, per the specificity finding above.

**Offline component accuracy** (from `read_cnn/pileup/`, the analysis the model was built on):

| Component | Metric | n |
|---|---|---|
| WGS hotspot scorer (`score_wgs.csv`) | ROC-AUC **0.82**; LIKELY-EDIT recall **0.92** @ precision **0.75** | 94 (53 pos / 41 neg) |
| Genome-wide PoN (`wgs_offtarget_worklist_pon.csv`) | LIKELY EDIT **276 → 185** after PoN; 118/185 on-target | 2737 candidates |
| ECS→WGS depth transfer (`check_transfer_ecs.csv`) | signal preserved at 30×: VAF median 0.242 → 0.236; 51/51 positives survive | 258 |

## Container

Every process runs in `ghcr.io/dhslab/docker-scge-offtarget:260710`, built from the
`docker-scge-offtarget/` folder in [dhslab-docker-images](https://github.com/dhslab/dhslab-docker-images).
It pins the versions the scripts and models need — notably scikit-learn 1.8.0, which is what
`wgs_shape_model.pkl` was trained under. If you ever change that pin, the build's own model-load
check (and the runtime guard in `bin/features.py`) will complain rather than quietly mis-score.
