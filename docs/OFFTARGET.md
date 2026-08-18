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
bnd_snapshots/         one PNG per JUNCTION in bnd_review_queue.tsv (not per row)
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

## Breakends

An indel is one cut healed badly. A **breakend** is two cut sites joined to each other — and on
this cohort that is overwhelmingly what survives triage: **23 of 25 queue rows are multi-cut
deletions**, two cuts from the same guide's own target set with the segment between them excised.
The pipeline has always emitted these. Until now it never labelled them and never drew them.

`bin/review_filter_bnd.py` applies the breakend analogue of the indel filter. Three rules, not
five, because a junction has two ends and no length spectrum:

**1. At least 3 supporting reads.** Breakend support is thin — the cohort median is 1 read — so
this is the single most discriminating cut available, and it does most of the work.

**2. Within 10 bp of a PAM position.** Same rule as indels — and the caller has already applied
the *same* 10 bp cutoff via `-d/--max-mutation-distance`, whose default is **10**, not 25;
`get_indels.nf` passes no override. The two tests are not quite identical (the caller keeps an
event if *either* end is within range, this rule tests one end), but on this cohort **no junction
has a cut distance above 10 at all**, so rule 2 can never fire. Treat "0 dropped" here as
structural, not as evidence the rule works.

**3. The breakpoint is not promiscuous.** A breakpoint partnering with many unrelated loci is an
alignment hub, not a junction. **On-target sites are exempt**, and that exemption is not optional:
a real Cas9 cut generates junctions to many places, so the true cut sites are among the most
promiscuous breakpoints in the cohort — 6 of the 10 breakpoints with ≥5 distinct partners are the
intended TRAC, TRBC1, TRBC2 and B2M sites. Applying rule 3 without the exemption deletes the real
edits.

There is no rule 4 here. The matched-control test is applied per event *upstream*: the caller's
`-x/--max-in-control` (default 0) drops any junction with control support before it ever reaches
`bnd_info`.

> ⚠️ Rule 3 is currently **non-binding** — on this cohort it drops nothing once rule 1 has run.
> It is kept because it costs nothing and is the only defence against an alignment hub that
> happens to be absent from the matched control. Do not read "0 dropped by rule 3" as broken.

### From 1,022 breakends to 8 junctions

The SV VCF holds 1,022 records and the queue holds 25 rows, and the gap between those numbers has
been asked about more than once. **It is one threshold, not a cascade.**

| stage | count | what removed the rest |
|---|---:|---|
| breakend entries in `bnd_info` (= VCF records) | **1,022** | — |
| after rule 1, `reads ≥ 3` | **25** | 997 junctions with 1–2 reads |
| after rules 2, 3, 4 | **25** | **nothing — all three dropped zero** |
| collapsed to events | **8** | grouping, not filtering |

Support is the whole story. The distribution of reads per junction:

```
reads:  1    2   3   4   5   6   7   8   9  13  15
count: 981  16   6   5   1   3   5   1   2   1   1
```

**981 singletons and 16 doubletons — those 997 are the entire reduction.** The VCF is the *ungated*
sibling of the queue: both read the same `bnd_info` column, and the VCF simply writes every entry
without applying `--min-reads`. So 1,022 vs 25 is not two analyses disagreeing; it is one analysis
before and after its only real threshold.

The last step, 25 → 8, is grouping and is described below — it is not a filter.

> **All 25 gated rows are `is_target == 1`**, and rules 3 and 4 exempt on-target sites. Both rules
> were therefore *structurally* inert on this cohort, not merely inactive. "0 dropped" is not
> evidence that they work.

### A queue row is not an event

**The 25 rows in `bnd_review_queue.tsv` are 8 junctions.** The caller reports each event from both
ends, at a few bp of position jitter, and under both strand orientations, so one deletion arrives
as three to five rows. Two consequences, and both have bitten:

- **Never quote a row count as an event count.** 25 breakends across 8 samples means 8 events, one
  per sample.
- **Never quote a row's `reads` as the event's support.** ARID4A's rows say 3–4 reads each; the
  junction carries 18 against 123–129× depth — a factor of 6.

`bin/bnd_snapshots.py` collapses rows to junctions on `(sample, sorted([bin, partner_bin]))` — bins
the filter already computes — and renders one figure each.

### Reading a junction figure

Each PNG in `review/bnd_snapshots/` is a to-scale schematic of the excision over a 2×2 grid: left
and right breakpoint across, edited sample over its matched unedited control down. Per-column
x-axes only — the two breakpoints have unrelated coordinates.

**The evidence is the green reads.** A read spanning a junction aligns in two pieces: a primary
clipped at one breakpoint and a supplementary segment at the partner, linked by an `SA` tag. Any
read whose `SA` lands near the partner locus is drawn green, so a junction reads as a stack of
green alignments all terminating on one base — present in the edited panel, absent from the
control directly below it. That pairing is what makes the figure self-adjudicating, exactly as in
the indel snapshots.

⚠️ **This evidence used to be invisible.** `bin/pileup_snapshot.py` dropped `is_supplementary`
unconditionally, which filtered out precisely the alignments that constitute a junction. It is now
opt-in per caller (`keep_supplementary`), so indel snapshots are unchanged. The counts involved are
not marginal: at ARID4A, 48 of 60 SA-tagged reads point at the partner locus; at PDCD4, 67 of 71.

Note the two read counts in a figure header are different statistics and will not agree: the
suptitle's *junction reads* is the caller's filtered support summed over the queue rows, while a
panel's *reads with SA at partner* is every SA-tagged read at that breakpoint. The first is the
conservative one.

**Not shown: a coverage drop.** Considered and rejected on two independent grounds — these events
are 6–20% VAF, so there is no visible dent, and the per-sample `tagged.bam` is region-limited, so
coverage between the two target islands is zero whether or not anything was deleted. The panel
would have read as a deletion that isn't there.

| parameter | default | what it does |
|---|---|---|
| `review_filter_bnd` | `true` | run the breakend filter at all |
| `review_bnd_min_reads` | `3` | rule 1 |
| `review_bnd_max_cut_dist` | `10` | rule 2 |
| `review_bnd_max_partners` | `5` | rule 3 |
| `review_bnd_bin_size` | `1000` | bp bin used to group breakpoints |
| `review_sv_noise` | `null` | external DRAGEN SV panel (BEDPE) |
| `review_bnd_snapshots` | `true` | render the junction figures |
| `review_bnd_snapshot_window` | `150` | bp either side of each breakpoint |
| `review_bnd_max_junctions` | `200` | safety cap on figures rendered |

⚠️ **Only `WGS_hg38_v3.1.0` is safe for `review_sv_noise`.** Measured: `IDPF_WGS v3.0.0` flags 25
of 25 real junctions and `FF_Heme v3.1.0` flags 24 of 25 — either erases the entire result.

## The noise model

Rule 4 asks whether a call is *statistically distinguishable from this locus's background*, not
whether it appears on a blacklist. That matters because a site can carry 1% background in the
control and 40% signal in the treated sample; an existence test throws that edit away.

    AQ = -10 log10 P(X >= k | n, background at this locus)

The background is an empirical-Bayes Beta posterior built from the control's own reads
(`bin/noise_model.py`). `review_noise_model = 'matched'` (the default) uses **the sample's own
matched control**, so it needs no cohort at all — which is the entire reason the old panel of
normals could be deleted. On the 32-sample CAR-T cohort the two are equal: **64/64 confirmed edits
retained, 9 human-rejected, precision 0.877** — that is the PoN-only arm (queue 88) against the
32×single-sample arm (queue 91), `docs/NOISE_MODEL_EXPERIMENT.md` arms 2 and 6. Equal, and the
matched model does it without a cohort.

> **Which precision figure to quote.** 0.877/9-rejected above is the *PoN-vs-matched equivalence*
> comparison. **Production is 0.889 with 8 rejected** (queue 96, 64/64 confirmed) — rule 1 on,
> matched AQ, caller-derived cut distance, `NOISE_MODEL_EXPERIMENT.md` "rule 1 ON (production)".
> Both numbers are measured on the same 32 tables and both are correct; they are different rows of
> the same experiment, and the retained-edit count (64/64) is identical in every arm. Always name
> the configuration — quoting either bare makes the two look like a contradiction, which is what an
> earlier version of this file and `README.md` between them managed to do.

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

## A second one, in the VCF output

Until now **every `*.offtarget_svs.vcf` this pipeline ever wrote contained zero records**.

`bin/bnd_from_indels_to_vcf.py` had its tab and newline escapes written doubled — `'\\t'` and
`'\\n'`, which in Python source is a backslash followed by a letter, not a control character. The
header line was split on the two-character string backslash-t, never matched, so the `bnd_count`
lookup failed, defaulted to `'0'`, and **every row took the `continue` branch**. The output was a
single line of literal backslash-n text; `wc -l` reported 0. Nothing downstream checked, because a
VCF with no records is exactly what a sample with no breakends should produce.

Nothing was miscalled and no analysis is affected — the triaged breakend queue is produced by a
different script from the same source column, and it was always correct. What was lost is a
published output file: the per-sample breakend VCF was empty for every sample ever run.

The fix also carries evidence that the old code discarded: `SR=` and `CTRL=` now report supporting
and control read counts, and reciprocal junctions are linked with `MATEID` rather than emitted four
times. On the 32-sample CAR-T cohort this is **1022 breakend records where there were 0**, verified
in the 2026-08-17 run: 32/32 samples populated, 3–95 records each.

Only the VCF step needs re-running, not the caller — the `bnd_info` column it reads from was never
wrong.

### What this fix does *not* do

It does not put breakends in the HTML report, and an earlier draft of this section wrongly said it
would. `COMPILE_REPORT_JSON` does consume the VCF, and the report JSON's `tables.bnd_vcf` is now
populated for all 32 samples — but **`bin/make_scge_report.qmd` never references `bnd_vcf`** (grep
it: zero hits), and neither does `bin/make_scge_excel.py`. The report's SV panel is driven by
`tables.on_target_sv_transgene`, a transgene-junction annotation table that is a different thing
entirely and is empty in **all 32** samples of this cohort — so that panel still reads "No on-target
SV/transgene data available in the report", and did so for its own unrelated reason all along.

Surfacing breakends in the report is therefore a template change, not a data change, and the data
is now sitting there waiting for it. Two separate defects; only the first is fixed.

> `bin/tsv_to_vcf.py` carries the same defect and is deliberately untouched: no module references
> it. Fix it before wiring it to anything.

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
