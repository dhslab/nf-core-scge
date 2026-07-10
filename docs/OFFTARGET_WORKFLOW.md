# Unified CRISPR Off-Target Workflow

A self-contained arm added to nf-core-scge on branch `feat/offtarget-wgs`. It runs
via a **named entry** (`-entry OFFTARGET`) and **does not touch the default SCGE pipeline**.

Per the PI (2026-07-09): the legacy `crispr_ml_*` read-level classifier is **deprecated** for
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

## What it produces

| Output (`<outdir>/offtarget/`) | From | Meaning |
|---|---|---|
| `wgs_offtarget_worklist_pon.csv` | WGS_WORKLIST → PON_OFFTARGET_FILTER | genome-wide, homology-free, PoN-filtered ranked off-target worklist |
| `offtarget_report.csv` | RECONCILE_OFFTARGET_REPORT | the worklist annotated with `is_hotspot` / `ecs_confirmed` / `ecs_if` |
| `<id>.offtarget_analysis.tsv` | ECS_INDELS | per-hotspot ECS error-corrected edit call + VAF (truth) |
| `training.tsv` | BUILD_TRAINING_TABLE | one row per (guide × hotspot): WGS features + ECS VAF + label (feeds the **offline** trainer) |
| `recall_vs_vaf.{csv,png}` | RECALL_VS_VAF | WGS recall of ECS-confirmed edits vs ECS VAF, and the trustworthy-VAF floor |

## Run modes (auto-detected from the samplesheet `datatype` column)

- **paired** (ecs + wgs rows): everything above — this is how the WGS-only model is trained/validated.
- **wgs_only**: worklist + PoN + report (no training/recall; empty ECS channels skip those steps).
- **ecs_only**: ECS truth tables only.

## Usage

```bash
nextflow run . -entry OFFTARGET -profile ris \
    --input offtarget_samplesheet.csv \
    --outdir ./results_offtarget
```
`run_offtarget.sh` wraps this in the RIS `bsub`. Samplesheet template:
`assets/offtarget_samplesheet_template.csv`.

**Samplesheet** (`sample,datatype,guide,edited_cram,control_cram,target_file,vcf`):
- `datatype` ∈ {ecs, wgs} — routes the row.
- `guide` — the join key between an ECS sample and its WGS counterpart.
- WGS `edited_cram` must be the DRAGEN `<base>_tumor.cram`; the matched normal `<base>.cram`
  and `<base>.hard-filtered.vcf.gz` must sit **beside it** (the scripts derive them by name).
- All paths **absolute** (CRAMs/VCFs are read from bind-mounted storage, not staged).

## Key params (in `nextflow.config`)

| param | default | purpose |
|---|---|---|
| `offtarget_shape_model` | `assets/models/wgs_shape_model.pkl` | pileup shape ranker (fixed asset; train offline) |
| `offtarget_min_af` | 0.05 | DRAGEN tumor AF floor for candidate indels |
| `offtarget_min_span` | 8 | spanning-read coverage gate |
| `offtarget_hi_score` | 0.60 | shape score ≥ this = "detected" for the recall curve |
| `offtarget_target_recall` | 0.80 | recall level whose VAF floor is reported |
| `offtarget_hotspot_pad` | 25 | bp window to match a worklist hit to a predicted hotspot |
| `offtarget_snapshots` | false | render IGV-style pileup PNGs for LIKELY EDITs |

## What this workflow is (and is not)

Scoped honestly against the real validation data (`read_cnn/pileup/`):

- **It is** a **hotspot edit-confirmation + genome-wide screen**. At known/nominated hotspots the
  WGS shape scorer recovers edits well (see numbers below); genome-wide it produces a PoN-filtered,
  ranked worklist for review.
- **It is not (yet)** a proven *homology-free de novo* off-target detector, and it does **not** yet
  carry a demonstrated detection floor below ~5% VAF. Treat WGS-only calls as trustworthy **≥5% VAF**;
  below that the workflow has no ground truth to stand on (see gaps).

## Measured accuracy (from `read_cnn/pileup/`)

| Component | Metric | n |
|---|---|---|
| WGS hotspot scorer (`score_wgs.csv`) | ROC-AUC **0.82**; LIKELY-EDIT recall **0.92** @ precision **0.75** | 94 (53 pos / 41 neg) |
| Genome-wide PoN (`wgs_offtarget_worklist_pon.csv`) | LIKELY EDIT **276 → 185** after PoN; 118/185 on-target | 2737 candidates |
| ECS→WGS depth transfer (`check_transfer_ecs.csv`) | signal preserved at 30×: VAF median full 0.242 → d30 0.236; 51/51 positives survive | 258 |
| Recall vs ECS VAF (`score_wgs.csv`) | 1.00 (0.05–0.20), 0.96 (0.20–0.50), 0.83 (>0.50) | 53 positives |

## Open gaps (do not overclaim past these)

1. **Sub-5% floor is unmeasured.** The truth set contains **zero** positives below VAF 0.05
   (min detected = 0.056), so the recall-vs-VAF curve cannot yet certify a floor under 5%. At 30×
   a 2% edit ≈ 0.6 supporting reads — physically near-unrecoverable — but that is *asserted*, not shown.
2. **Off-target discovery has no positive control.** The genome-wide arm finds no convincing novel
   off-targets; the residual LIKELY-EDIT hits recur at TCR loci (chr14:22.5M / chr7:142.8M) and are
   lineage/mapping artifacts, not guide off-targets. "Finds nothing" is validated; "would fire on a
   real off-target" is not.
3. **The paired Nextflow glue arm** (`hotspot_to_table.py`, `join_training_table.py`,
   `recall_vs_vaf.py`) is validated only at the script/schema level — the first RIS run
   (`run_offtarget_aavs1.sh`, AAVS1 paired subset) is what proves the (guide, chrom, start) join
   end-to-end (an ECS/WGS coordinate off-by-one yields an empty `training.tsv`; the joiner warns).

## Container note

Processes use `ghcr.io/dhslab/docker-scge:latest`, which must contain the new `bin/` scripts and
their deps (pysam, scikit-learn matching `wgs_shape_model.pkl`, joblib, pandas, numpy, **matplotlib**).
`docker-scge` already carries the ML stack for the existing pipeline; confirm `matplotlib` is present
(needed for snapshots + the recall curve) and rebuild the image so the new scripts are baked in
(the scripts cross-import, so they must be co-installed on `PYTHONPATH`, as the pipeline already does
for `crispr_ml_features.py`).
