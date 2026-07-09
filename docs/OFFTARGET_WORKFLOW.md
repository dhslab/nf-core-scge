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

## Confidence / validation status

- **WGS worklist + PoN** (`worklist_from_vcf.py`, `pon_filter.py`): the validated core. Established
  result to reproduce: LIKELY EDIT 276→185 after PoN; PLCB2 positive control survives; zero
  convincing homology-blind off-targets.
- **Paired training/recall arm** (`hotspot_to_table.py`, `join_training_table.py`, `recall_vs_vaf.py`):
  **new code, not yet run on real data.** Python compiles and the schema-level logic is correct
  against a real `offtarget_analysis.tsv`, but validate on the first RIS run — in particular the
  (guide, chrom, start) join (an ECS/WGS coordinate off-by-one would yield an empty `training.tsv`;
  `join_training_table.py` warns if so).

## Container note

Processes use `ghcr.io/dhslab/docker-scge:latest`, which must contain the new `bin/` scripts and
their deps (pysam, scikit-learn matching `wgs_shape_model.pkl`, joblib, pandas, numpy, **matplotlib**).
`docker-scge` already carries the ML stack for the existing pipeline; confirm `matplotlib` is present
(needed for snapshots + the recall curve) and rebuild the image so the new scripts are baked in
(the scripts cross-import, so they must be co-installed on `PYTHONPATH`, as the pipeline already does
for `crispr_ml_features.py`).
