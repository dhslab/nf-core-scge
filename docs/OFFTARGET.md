# CRISPR off-target detection

Everything about finding off-target edits lives here. There are **two separate things** in the
pipeline, and mixing them up is the most common confusion:

| | what it is | when you use it | how to run |
|---|---|---|---|
| **Review filter** | turns a routine run's call table into a short list a human can actually read | every normal run | on by default |
| **`-entry OFFTARGET`** | a two-assay (ECS + WGS) investigation arm | building or checking the model | `-entry OFFTARGET` |

Most people only ever need the first one. Start there.

New to this? [`OFFTARGET_DEMO.md`](OFFTARGET_DEMO.md) runs the whole thing on a 5.5 KB toy genome
in about a minute — no cluster, no cohort data.

---

# Part 1 — Automated review (the routine path)

The pipeline predicts thousands of candidate sites per guide and scores them all. It used to hand
back a list you had to open one by one in IGV. Now five rules cut that down and draw you a picture
of each survivor.

**On the 25-sample CAR-T WGS cohort: 238 sites to review became 62, keeping all 61 real edits.**

| | queue | real edits | false positives | recall | precision |
|---|---|---|---|---|---|
| before | 238 | 61 | 177 | 1.000 | 0.256 |
| now | **62** | **61** | **1** | **1.000** | **0.984** |

## What you get

Run the pipeline normally, then look in `<outdir>/review/`:

```
review_queue.tsv       the sites to actually look at
review_queue_all.tsv   every gated site + why_dropped  (the audit trail)
bnd_review_queue.tsv   the same triage applied to breakends
snapshots/             one PNG per site in review_queue.tsv
```

Each snapshot has two panels: **the edited sample on top, its matched unedited control below** —
same locus, same scale. That pairing is the whole point. Germline variants and alignment artifacts
show up in *both* panels; a real edit shows up in only one. Most calls are settled by glancing at
the bottom panel.

Nothing is silently thrown away. `review_queue_all.tsv` lists every site that cleared the gate
along with the reason it was dropped, so you can check any decision without re-running anything.

## The five rules

Sites enter the queue at `indel_reads >= 10 AND indel_fraction >= 0.05`. Then:

**1. The matched control must be clean** (control VAF < 5%)
If the variant is in the unedited sample from the same donor, Cas9 didn't make it. Biggest filter
by a wide margin.

**2. The indel sits within 10 bp of a PAM**
Cas9 cuts ~3 bp from the PAM. Confirmed edits sit at a median of 2–4 bp; background indels don't.
Depth-stable — resample a real edit down to 3 reads and 100% still pass.

**3. At least 3 distinct indel lengths**
NHEJ makes a spectrum of deletion sizes at one cut site. One length repeated across every read is
an alignment artifact.
⚠️ **Depth-sensitive.** At 3 supporting reads only ~49% of real edits still pass. Safe here because
everything in the queue already has ≥10 reads. Don't reuse this rule to chase lower-VAF events — it
quietly reimposes a read-count floor.

**4. Not background noise** — see [the noise model](#the-noise-model).

**5. Not in a repeat region**
The candidate panel is built from sequence homology, so it's deliberately enriched for repetitive
and paralogous sequence — exactly where aligners invent indels.

On-target sites are exempt from rules 4 and 5, since they're shared across samples by design.

## The noise model

Rule 4 asks whether a call is *statistically distinguishable from this locus's background*, not
whether it appears on a blacklist. That matters because a site can carry 1% background in the
control and 40% signal in the treated sample; an existence test throws that edit away.

    AQ = -10 log10 P(X >= k | n, background at this locus)

The background is an empirical-Bayes Beta posterior built from the control's own reads
(`bin/noise_model.py`). `review_noise_model = 'matched'` (the default) uses **the sample's own
matched control**, so it needs no cohort at all — which is the entire reason the old panel of
normals could be deleted. On the 32-sample CAR-T cohort the two are equal: 64/64 confirmed edits
retained, 9 human-rejected, precision 0.877, and the model does it without a cohort.

`review_aq_min` (default 5) is the cut. Measured with the depth floor on, 3-8 all reproduce the
PoN result exactly; 10 starts costing confirmed edits.

**The depth floor.** A control with no alt reads at depth *d* shows the background is below ~1/*d*,
not that it is zero. Left alone the posterior collapses to 3.5e-5 against a 170x control, about
100x beyond what 170 reads can support, and every call at a locus the control never sampled deeply
starts looking significant. `review_depth_floor` (default true) bounds the posterior mean at
1/control_depth. It is harmless at the current `VAF >= 0.005` gate and essential if that gate is
ever lowered for a high-sensitivity run.

**Germline falls out for free.** A germline het sits near 50% in the matched control, so a 50%
observation in the treated sample is unsurprising and scores low. Measured: 0 of 132 germline-like
sites survive the cut, with no dedicated germline rule.

**Single guide and no noise model?** Rule 4 falls back to cross-guide recurrence, which needs >=2
guides passed in ONE invocation and therefore **cannot fire for a single guide**. The filter prints
a boxed `RULE 4 IS NOT ACTIVE` warning in that state; `--strict-fallback` /
`review_strict_fallback` makes it a hard failure instead.

**No unedited control at all?** Use rule 5 alone (`--review_repeat_beds`). Repeat annotation needs
no controls, no cohort and no guide context, and `--snv-noise` adds the external DRAGEN panel.
Both are weaker than the noise model — the DRAGEN panel has an indel-capable record for only 3.6%
of queried loci, so an arbitrary floor decides the rest — but they work on day one.

**What about the GATK panel of normals?** Tested and rejected. `1000g_pon.hg38.vcf.gz` catches 15
of our artifacts but also hits **4 of the 61 real edits** — it's a Mutect2 SNV panel from blood
normals, not an indel-artifact map, so it costs recall:

| resource | artifacts caught (of 177) | real edits wrongly hit (of 61) |
|---|---|---|
| GATK 1000g PoN | 15 | **4** |
| RepeatMasker | 100 | 0 |
| Tandem repeats (TRF) | 54 | 0 |
| this run's own PoN | 163 | 0 |

## Is there a model here?

There's a trained one, and it is **not** what runs. A gradient-boosted classifier on the same
features reaches out-of-fold AUC 0.9992 under leave-one-sample-out validation, and at 100% recall
returns the same 62 sites with 1 false positive that the rules reach. It wins nothing worth a
pickle file and a scikit-learn version pin, so the rules ship.

One finding from building it is worth keeping: **never give a model `is_target`.** The truth set is
52 confirmed on-target edits against 1 confirmed off-target, so the model learns "on-target ⇒ real"
and scored the cohort's only genuine off-target at p = 0.113 — it nearly threw away the one event
the assay exists to find. Dropping that feature raised AUC *and* moved that event to p = 1.000.

## Parameters

| parameter | default | what it does |
|---|---|---|
| `review_filter` | `true` | run the filter at all |
| `review_min_reads` / `review_min_vaf` | `2` / `0.005` | the entry gate |
| `review_max_control_vaf` | `0.05` | rule 1 |
| `review_max_cut_dist` | `10` | rule 2 |
| `review_min_distinct_len` | `3` | rule 3 |
| `review_noise_model` | `matched` | rule 4: beta-binomial vs the sample's own control |
| `review_aq_min` | `5` | AQ below this is background |
| `review_depth_floor` | `true` | bound the posterior at 1/control_depth |
| `review_strict_fallback` | `false` | fail, not warn, when rule 4 has no source |
| `review_snv_noise` | `null` | external DRAGEN noise panel for rule 6 |
| `review_repeat_beds` | `null` | comma-separated repeat BEDs for rule 5 |
| `review_snapshots` | `true` | render the snapshot packet |

Run it standalone on existing call tables:

```bash
bin/review_filter.py results/*.offtarget_analysis.tsv \
    --noise-model matched \
    --repeats /path/rmsk.no_simple.bed /path/trf.bed \
    --min-reads 2 --min-vaf 0.005 \
    -o review_queue.tsv
```

Add `--keep-all` to get every gated row with its `why_dropped` reason instead of a filtered list.

## Two caveats before quoting a number

**Thresholds were fitted on one cohort.** 10 bp and 3 distinct lengths come from 25 CAR-T WGS
samples, not from first principles. They hold across every sample in that cohort, but a new guide
panel deserves a re-check against its own controls.

**This only operates above the existing gate.** It automates *review* of sites clearing
`reads >= 10 & VAF >= 0.05`; it does not change what's *detected*. At the cohort's ~207× depth that
gate is a genuine ~5% VAF detector. Going below 5% is a different problem this filter doesn't touch.

## A bug worth knowing about

Until recently `find_edited_reads.py` reported **no control support at any site, in every run ever
done**. `-x/--max-in-control` drops events seen in the matched control — correct, and still how it
works — but the per-site summary was totalled *after* that filter had already removed the very
events carrying the evidence, so `control_indel_reads` was always 0.

Nothing was miscalled. But rule 1, the most powerful filter here, was a no-op against those tables.
The fix moves the summation before the filter; old and new output differ in exactly two columns
(`control_indel_reads`, `control_indel_fraction`) and are byte-identical everywhere else.

**Call tables made before the fix still have zeroed control columns.** `review_filter.py` detects
this and warns. Re-run the caller to get rule 1 back.

---

# Part 2 — `-entry OFFTARGET` (the investigation arm)

A separate arm that runs via a named entry and does not touch the default pipeline. It pairs deep
error-corrected sequencing (**ECS — the truth arm**) against ordinary WGS (**the arm under test**)
to ask: could WGS alone have found these edits?

## Running it

```bash
sbatch run_offtarget.sh --input <samplesheet.csv> --outdir <dir> [--snapshots]
```

or directly (from somewhere that can `sbatch` — not the interactive exec node):

```bash
nextflow run . -entry OFFTARGET -profile ris2,apptainer \
    --input offtarget_samplesheet.csv --outdir ./results_offtarget -resume
```

On Compute1 (LSF) use `-profile ris` under `bsub`.

## Samplesheet

Columns: `sample,datatype,guide,edited_cram,control_cram,target_file,vcf`
(template at `assets/offtarget_samplesheet_template.csv`)

- `datatype` is `ecs` or `wgs` — this sorts each row into the right arm.
- `guide` links an ECS sample to its WGS partner; same guide means same experiment.
- For a WGS row, point `edited_cram` at the DRAGEN `<name>_tumor.cram` and keep `<name>.cram` (the
  normal) and `<name>.hard-filtered.vcf.gz` beside it. The scripts find them by name.
- Use absolute paths — CRAMs are read off storage, not copied in.

What runs depends on what you give it: **ECS + WGS** → everything; **WGS only** → worklist and
report; **ECS only** → hotspot truth tables. The run logs which mode it picked, and fails
immediately with the offending file and row if a DRAGEN sidecar is missing.

## What you get

Everything lands in `<outdir>/offtarget/`:

- **`wgs_offtarget_worklist_pon.csv`** — the main result. Every candidate edit found in the WGS,
  ranked, germline and artifacts already knocked out by the PoN. Start here.
- **`offtarget_report.csv`** — the worklist plus whether each hit is a known hotspot and whether
  ECS confirmed it.
- **`<id>.offtarget_analysis.tsv`** — the ECS answer at each hotspot: edited or not, at what VAF.
- **`recall_vs_vaf.csv` / `.png`** — how often WGS catches an ECS-confirmed edit, by VAF.
- **`offtarget_metrics.json` / `.txt`** — PR-AUC and F2/F5 for the shape score. **Read the
  denominator section first** — precision here is against the ECS label, *not* human review.
- **`training.tsv`** — WGS features lined up against the ECS answer. Only used for retraining.

The last three need both ECS and WGS.

## Optional review aids

**Snapshots** (`--offtarget_snapshots true`) — an IGV-style pileup per LIKELY EDIT, edited beside
matched normal. Off by default since it renders one image per hit.

**Per-read tags** (`--offtarget_tagged_bam true`) — writes `<id>.tagged.bam` where every read
carries an `XC` tag naming how the caller classified it. In IGV: load the BAM → right-click →
**Color alignments by → tag → `XC`**.

| tag | meaning |
|---|---|
| `Edited_Deletion_<N>bp` / `_Insertion_` / `_Duplication_` | the event that was called |
| `Edited_BND_<chrom>` | breakend — the read's mate maps to `<chrom>` (e.g. the transgene) |
| `Edited_SoftClip` | recovered by realigning a soft clip |
| `Unedited_WT` | spans the target cleanly — this is the denominator |
| `Skipped_*` | excluded by a read filter, shown so you can see *why* |

Each alignment appears exactly once, so IGV's depth is real.

**Don't expect tag counts to match the TSV.** Tags are per *alignment record*; `indel_reads` is per
*fragment* (R1/R2 collapsed) and taken after site-level filters — so `Edited_*` records run roughly
2× higher where both mates cover the cut.

⚠️ **Keep it off for cohort sweeps.** Measured on one AAVS1 ECS sample (1,149 targets, ~11,000×):
800 MB BAM, peak RSS 1.4 GB → 6.3 GB, wall clock ~30 → ~46 min. The process gets 24 GB instead of
8 GB when the flag is set. Check one sample before enabling across a cohort.

## The shape model

`assets/models/wgs_shape_model.pkl` — a `HistGradientBoostingClassifier` (`max_depth=3`), not a
neural network. Role: `stage2_shape_ranker`. Trained on 51 positive / 383 negative loci from the
CAR-T ECS cohort, depth-augmented at 0/20/30/50× so it behaves at 30× WGS.

**The deployed model uses 7 features:**
`indel_frac`, `conc_ratio`, `pos_conc`, `pos_mad`, `modal_len`, `modal_mapq`, `softclip_frac`.

> ⚠️ **`features.py:MODEL_FEATURES` lists 10**, adding `cut_dist`, `homopolymer_len` and
> `repeat_frac`. Those three were built and evaluated but the 10-feature model was **never
> shipped** — it bought 2.1× precision at the cost of CAR-T recall (52/52 → 47/52). `score.py`
> selects inputs by `bundle["features"]`, so the deployed pickle stays at 7 regardless. **If you
> retrain today you will get a 10-feature model**, which is a different model from the one these
> numbers describe. The ideas weren't wasted: `cut_dist` and repeat context are now rules 2 and 5
> of the review filter, where they're measurable and explainable.

Two weaknesses worth volunteering before someone else finds them:

- **51 positives is small.** Treat it as a prioritiser, not an oracle.
- **It's a ranker, not a detector** — it demotes some unambiguous edits. Which is why the rescue
  exists.

**The high-evidence rescue.** Scoring on the model alone recovered 47/52 confirmed CAR-T edits. The
five misses weren't close calls — indel fractions 0.25–1.00, clean controls, 130–244 spanning reads,
pileups a reviewer calls instantly — that the ranker scored 0.27–0.50. So `score.py` calls LIKELY
EDIT on unambiguous evidence regardless of model score, via two arms, both gated *behind* the
low-MAPQ and matched-normal vetoes. The rescue can override the **model**, never the **evidence**,
so it can't resurrect a germline variant:

1. **under-scored** — `indel_frac ≥ 0.15`, `conc_ratio ≥ 0.5`, `spanning ≥ 20`
2. **multi-allelic** — model score ≥ `offtarget_hi_score` with the same depth/burden evidence,
   waiving the clonality requirement

`call_basis` records which fired. Across 99,308 scored CAR-T rows this promoted exactly 6 sites —
all 6 the confirmed edits, taking recall to **52/52**. Set `--offtarget_rescue false` for
model-only behaviour.

## Retraining (`-entry TRAIN`)

The deployed model is a fixed asset; training is deliberately not in the OFFTARGET DAG.

```bash
# 1. paired run -> labeled training table
sbatch run_offtarget.sh --input paired_samplesheet.csv --outdir results
# 2. fit
nextflow run . -entry TRAIN -profile ris2,apptainer \
    --input results/offtarget/training.tsv --outdir results
# 3. deploy
sbatch run_offtarget.sh --input cohort.csv --outdir results2 \
    -- --offtarget_shape_model results/train/wgs_shape_model.pkl
```

It runs in the off-target container so the pickle is written under the same scikit-learn the
pipeline scores with. See the 7-vs-10 feature warning above before you deploy the result.

## Reading the metrics (the denominator matters more than the number)

There are **two denominators**, kept strictly apart:

| | denominator | what it supports | where |
|---|---|---|---|
| 1 | ECS label (a real two-class label) | PR-AUC, precision, F-beta | `offtarget_metrics.{json,txt}` |
| 2 | human manual review (positives only) | **recall only** — the gold standard | `bin/validate_recall.py` |

**You cannot compute precision from manual review.** It has 55 confirmed positives and 80,440 NaNs,
and NaN means *not reviewed*, not *reviewed and rejected*. There are no confirmed negatives in it.
Filling those NaNs with 0 would count every genuine discovery the reviewers never reached as a
false positive. That's why `validate_recall.py` reports recall only, and stays that way.

**Why recall reads low against the ECS label.** The raw ECS label counts any nonzero indel fraction
as an edit, and at ECS depth that's overwhelmingly noise — most "edits" sit below 0.5% VAF with a
median of 3 indel reads out of ~5,000. Of the high-VAF ones, most carry the same indel in the
matched normal, i.e. they're **germline**. The WGS arm rejects those correctly. So the credible
denominator requires somatic status, a VAF floor and a read-count floor; sites with no spanning WGS
reads are reported as **unevaluable, not missed**.

Report PR-AUC and recall-weighted F-beta, not ROC-AUC — at this prevalence ROC-AUC is dominated by
true negatives and reads flatteringly high. F-beta is monotone in beta, so F1 is always an
*endpoint* of {F1, F2, F5}. **F1 in the middle means the weights got applied backwards.**

## What this arm is, and isn't

- **It is** a hotspot edit-confirmation plus genome-wide screen. On-target recovery is proven
  end-to-end: on the first real AAVS1 run the on-target came back from WGS alone.
- **It also** recovers the one confirmed *off*-target we have — PLCB2 chr12:32,679,410 (90% VAF,
  ECS-confirmed) — the proof WGS-only can find a real homology-based off-target.
- **It is not** a workflow with a quantified low-VAF detection floor. Across both cohorts the
  human-reviewed truth set holds only **2 confirmed off-targets, both >85% VAF**. There's no
  low-abundance off-target population to measure against. That's editing being highly
  on-target-specific, not a bug — but it means sub-5% sensitivity is unproven. **Trust WGS-only
  calls at hotspots ≥5% VAF.**

---

## Known gaps, stated once

- **No curated low-VAF positives**, so the model's core claim is untested.
- **Sub-5% VAF at 30× WGS is physics**, not tuning — a 1% VAF edit yields ~0.3 supporting reads.
- **No de-novo off-target positive control**: the genome-wide LIKELY EDITs that aren't predicted
  hotspots are unreviewed. (Of 49 genome-wide calls, 9 fall on predicted hotspots and all 9 are
  ECS-confirmed.)
- **Base-editor (CBE/ABE) substitutions are not detected.** This pipeline calls indels only.

## Container

Every process runs in `ghcr.io/dhslab/docker-scge-offtarget:260710`, built from
`docker-scge-offtarget/` in [dhslab-docker-images](https://github.com/dhslab/dhslab-docker-images).
It pins scikit-learn 1.8.0, which is what `wgs_shape_model.pkl` was trained under. Change that pin
and the build's own model-load check (plus the runtime guard in `bin/features.py`) will complain
rather than quietly mis-score.
