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
- **`training.tsv`** — one row per hotspot lining up the WGS features against the ECS answer. Only
  used to retrain the model offline; ignore it on a normal run.

The last two only show up when you give it both ECS and WGS.

## What runs depends on the samplesheet

The `datatype` column decides:

- **ECS + WGS rows** → the whole thing: worklist, report, training table, recall curve. Use this to
  build or check the model.
- **WGS only** → worklist and report. No training table or recall curve.
- **ECS only** → just the hotspot truth tables.

## Running it

```bash
nextflow run . -entry OFFTARGET -profile ris \
    --input offtarget_samplesheet.csv \
    --outdir ./results_offtarget
```

`run_offtarget.sh` does the same thing under `bsub` on RIS, and there's a filled-in example at
`assets/offtarget_samplesheet_template.csv`.

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

## Container

Every process runs in `ghcr.io/dhslab/docker-scge-offtarget:260710`, built from the
`docker-scge-offtarget/` folder in [dhslab-docker-images](https://github.com/dhslab/dhslab-docker-images).
It pins the versions the scripts and models need — notably scikit-learn 1.8.0, which is what
`wgs_shape_model.pkl` was trained under. If you ever change that pin, the build's own model-load
check (and the runtime guard in `bin/features.py`) will complain rather than quietly mis-score.
