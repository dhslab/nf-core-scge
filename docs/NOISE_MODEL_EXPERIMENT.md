# Can a probabilistic noise model replace the panel of normals?

Measured 2026-08-13 on the 32-sample CAR-T WGS cohort (`results_cart_ponfix`, 99,238 rows,
1,498 clearing the `reads>=2, VAF>=0.005` gate).

**Short answer: yes, but not with the external DRAGEN panel.** A beta-binomial test against the
sample's *own matched control* reproduces the PoN's discriminative power exactly, with no cohort
required. The DRAGEN systematic-noise panel does not come close, and for a reason that is
structural rather than tunable.

## Truth and its limits

Two sources, kept separate because they answer different questions.

- **Curated WGS label** — `Manual Indel Review/cart_wgs/cart_wgs_merged.xlsx`. Inside the stratum
  `indel_fraction >= 0.05 AND indel_reads >= 10` every row was adjudicated by eye, so a blank
  `manual_review` there is a **rejection**, not an absence. 241 rows → 41 confirmed / 138 rejected
  truth keys. Two-class, but only at VAF ≥ 5%.
- **ECS review** — `cart_ecs/cart_ecs_merged.csv.gz`, 56 confirmed edits at ~1,960× depth,
  reaching down to 2.48% VAF. **Positives only** (NaN = not reviewed), so it supports recall and
  nothing else — but it is the only truth that reaches below the WGS review's 5% floor.

Two caveats that affect every number below:

- Truth is keyed by `(guide, chrom, start)`, so one truth key joins to several scored rows when a
  guide has replicate samples. That is why 41 confirmed keys appear as 64 confirmed scored rows.
- 5 of 179 truth keys carry conflicting labels across replicates (2.8%); the key takes the max.

## Part A — the PoN-off experiment

Two traps, both real: the script's own defaults are `10 reads / 5% VAF` while the pipeline passes
`2 / 0.005`, and **omitting `--pon` does not disable rule 4** — it silently falls back to
cross-guide recurrence. Rule 4 is only truly off with `--max-guides 9999`.

| # | Configuration | Queue | conf | rej | precision |
|---|---|---:|---:|---:|---:|
| 1 | PoN + DRAGEN panel (**what ships today**) | 87 | 64 | 8 | 0.889 |
| 2 | PoN only | 88 | 64 | 9 | 0.877 |
| 3 | cross-guide fallback + panel | 88 | 64 | 9 | 0.877 |
| 4 | rule 4 OFF + panel (**rule 6 standing alone**) | 91 | 64 | 9 | 0.877 |
| 5 | rule 4 OFF, no panel (the floor) | 92 | 64 | 10 | 0.865 |
| 6 | 32 × single-sample, no PoN (**the target user**) | 91 | 64 | 9 | 0.877 |

Read off the table:

- **Rule 6 standing alone removes exactly one site** (92 → 91) — the same site it removes when the
  PoN is present (88 → 87). The PoN removes four *different* sites (92 → 88). The two are
  **disjoint**: the external panel does not find what the PoN finds.
- The PoN's raw attribution says 10 drops, but disabling it only costs 4 queue rows, because the
  repeat rule reclaims 6 of them (rule 5 drops rise 9 → 15). First-match-wins attribution
  overstates any single rule's unique contribution.
- **Arm 6 = arm 4 exactly (91).** Rule 4 is the only cohort-dependent rule, so a single-sample
  submission loses precisely the PoN and nothing else.

## Part B — the beta-binomial scorer

`bin/noise_model.py`. `AQ = -10 log10 P(X >= k | n, background)`, background from an
empirical-Bayes Beta prior fitted over every control observation (here Beta(0.0043, 1.774), mean
0.244%) and updated per locus. An unobserved locus keeps the prior — which is what makes the
"no data here" case principled instead of `p = 0`.

Scored on the 232 rows carrying a curated label:

| baseline | AUC | AP | AQ≥30 recall | AQ≥30 precision |
|---|---:|---:|---:|---:|
| **matched** (own control) | **0.977** | **0.944** | 64/64 | 0.681 |
| both | 0.971 | 0.908 | 64/64 | 0.674 |
| loo (cohort controls) | 0.936 | 0.801 | 51/64 | 0.739 |
| panel, `p = MAX` | 0.946 | 0.872 | 64/64 | 0.312 |
| panel, `p = MEAN × N/NR` | 0.945 | 0.869 | 64/64 | 0.312 |
| panel, `p = MEAN` | 0.930 | 0.847 | 64/64 | 0.306 |

### Why the panel arms fail

**Only 54 of 1,498 gated loci (3.6%) have any indel-capable panel record.** The other 96.4% are
scored against `--floor-p`, so the floor *is* the filter. Sweeping it proves the point:

| `--floor-p` | AUC | AQ≥30 recall | precision |
|---|---:|---:|---:|
| 0.0001 | 0.954 | 64/64 | 0.283 |
| 0.001 | 0.945 | 64/64 | 0.312 |
| 0.01 | 0.928 | 63/64 | 0.423 |
| 0.05 | 0.918 | 56/64 | 0.629 |

That swing dwarfs the choice of `--panel-p` (ΔAUC 0.016 across all three). The `MEAN`-vs-`MAX`
distinction is real — `MEAN` is diluted by the full panel, verified at MAX/MEAN = 45.99 for NR=1
records, exactly the 46-sample panel size — and `MEAN` is indeed the worst of the three. But it is
a second-order effect when the panel only speaks for 3.6% of sites.

### Sub-5% recall, the thing most at risk

12 ECS-confirmed edits below 5% VAF joined the scored set. At AQ ≥ 30 the **matched** baseline keeps
**12/12**; the cohort baseline keeps 10/12. So the conservative test does *not* destroy low-VAF real
edits — the concern that motivated this check does not materialise for the matched baseline.

### Germline

Germline-like sites (matched-control VAF ≥ 30%), n = 132:

- matched baseline: **0/132** reach AQ ≥ 30 — germline is filtered with no dedicated rule, exactly
  as the binomial argument predicts.
- cohort baseline: **1/132** escapes. The cross-donor failure mode is real but rare here, because
  pooling 31 donors' controls does capture most germline. Worked example —
  `RXRB-KO-DNA chrX:105,625,831`, treated 9/14 (0.64 VAF), matched control 28/64 (0.44):
  AQ 9.7 matched vs **31.2** cohort.

## The decisive result

Applying the scorer *in place of rule 4* — arm 5 (rules with rule 4 disabled) plus an AQ cut:

| | Queue | conf | rej | precision |
|---|---:|---:|---:|---:|
| arm 5, no filter | 92 | 64 | 10 | 0.865 |
| **arm 5 + matched AQ ≥ 30** | **91** | **64** | **9** | **0.877** |
| arm 2 (the PoN), for reference | 88 | 64 | 9 | 0.877 |
| arm 5 + matched AQ ≥ 60 | 77 | 59 | 6 | 0.908 |
| arm 5 + matched AQ ≥ 100 | 59 | 47 | 1 | 0.979 |

**AQ ≥ 30 on the matched-control baseline matches the panel of normals exactly** — same recall
(64/64), same rejected count (9), same precision (0.877) — while needing no cohort at all. Pushing
to AQ ≥ 60 buys precision 0.908 but costs 5 confirmed edits; that trade is not worth taking.

## What this means for the single-sample user

The distinction that matters is **matched control ≠ panel of normals**:

- A user with **one edited sample and its paired unedited normal** can drop rule 4 entirely and lose
  nothing. This is the realistic case and the recommendation.
- A user with **no control material whatsoever** falls back to the DRAGEN panel, where precision is
  0.31 rather than 0.68 and the result is governed by an arbitrary floor. The external panel is a
  weak backstop, not a PoN substitute.

## The failure mode Tier 1 actually has (and it is not missing controls)

All 32 samples carry a matched control: median depth **170×** (range 99–172), and none shows the
degenerate "control depth present but zero alt reads anywhere" signature left by the pre-fix caller.
So the no-control tier does not arise for this assay.

The real risk is the opposite of over-conservatism — it is **over-confidence on a clean control**.
66.4% of gated rows (995/1,498) have *zero* alt reads in the matched control. For those the
posterior collapses to a background of **3.5e-5**, roughly 70× below the fitted global prior mean of
2.44e-3. A clean control at 170× can honestly support a bound of about 1/170 = 0.59%; the posterior
instead asserts something near 1-in-28,000, because `alpha` stays at `a0` while `beta` grows with
depth.

On this cohort that never bites, and the reason is the gate, not the model: only **2** gated rows sit
below the control's 0.59% resolution limit, and **neither reaches AQ ≥ 30**. The `VAF >= 0.005` gate
is what protects the result.

That makes Tier 1 safe as currently configured and **fragile if the gate is ever lowered** — which is
exactly what a deliberate low-VAF off-target hunt would do. The fix is to floor the posterior
background at roughly `1 / control_depth` rather than let it run to `a0 / (a0 + depth)`, so the model
cannot claim more resolution than the control has. Not implemented here.

## What was implemented

The floor is on by default, and it **rescales AQ** — it does not cost anything at the operating
point, but it moves the threshold. Re-measured with the floor active:

| | Queue | conf | rej | precision |
|---|---:|---:|---:|---:|
| rule 4 disabled, no AQ | 92 | 64 | 10 | 0.865 |
| + floored AQ ≥ 3 | 91 | 64 | 9 | 0.877 |
| **+ floored AQ ≥ 5** (default) | **91** | **64** | **9** | **0.877** |
| + floored AQ ≥ 8 | 90 | 64 | 9 | 0.877 |
| + floored AQ ≥ 10 | 87 | 63 | 9 | 0.875 |
| the PoN, for reference | 88 | 64 | 9 | 0.877 |

3–8 form a plateau that reproduces the PoN exactly; 10 is the cliff where a confirmed edit is lost.
Default `review_aq_min = 5` sits in the middle. (Without the floor the equivalent threshold was 30 —
the same decision, different scale.)

Shipped as:

- `bin/noise_model.py` — `control_posterior`, `apply_depth_floor`, `fit_global_prior`, reused by
  the filter rather than duplicated.
- `bin/review_filter.py` — `--noise-model {off,matched,loo,both}`, `--aq-min`, `--no-depth-floor`,
  `--strict-fallback`. When the model is on it **replaces** rule 4 and the PoN is not passed
  alongside it, so a run's provenance stays readable. `why_dropped` reads
  `indistinguishable from control noise (AQ<5)`.
- Pipeline params `review_noise_model`, `review_aq_min`, `review_depth_floor`,
  `review_strict_fallback`, wired through `modules/local/review_filter.nf` into **both**
  invocations and registered in `nextflow_schema.json`.

**Default is `off`.** Verified: with defaults the 88 and 87 baselines still reproduce, and the
87-row queue is byte-identical to the pre-change output.

### The Tier 2 guardrail

The dangerous state is silent: no PoN, no noise model, and a single guide in the invocation means
cross-guide recurrence cannot fire, and the unfiltered queue is indistinguishable from a precise
one. That now prints a boxed `RULE 4 IS NOT ACTIVE` warning naming the three fixes, and
`--strict-fallback` / `review_strict_fallback = true` turns it into a non-zero exit. It stays silent
when `--noise-model` covers rule 4.

## Incidental, and worth shipping

Restricting the 1 GB panel to the positions a cohort can query collapses it from **97,851,753 to
4,891 records (60 KB)** and produces **byte-identical** filter output. That removes the ~6 min
streaming cost and makes the tabix question moot.

## Annotate in the caller, filter in the TSV

With the PoN gone, nothing in the filter needs a cohort — so the two-stage split had to be
re-justified. It survives, on cost asymmetry:

| stage | cost |
|---|---|
| `REVIEW_FILTER`, whole cohort, from TSVs | **1.78 s** |
| `GET_INDELS`, per sample, from CRAM | **50 min – 1h 28m** |

~2,000×. Every threshold here (`min_reads`, `min_vaf`, `max_cut_dist`, `min_distinct_len`,
`aq_min`) was calibrated by re-running the filter; moving those decisions into the caller converts
each experiment from seconds into hours and deletes the `why_dropped` audit trail — including the
**753 of 1,498** rows sitting 11–25 bp from the cut, which are exactly the evidence for the 10 bp
threshold. The intermediate TSVs are **27 MB** total, so there is no I/O argument either way.

So: computation moves to the caller, decisions stay in the filter. `add_features` now **consumes**
caller-supplied columns and derives them only as a fallback, which keeps tables from older callers
working unchanged.

### The cut-distance trap

`min_cut_distance` (caller) and `cut_dist_min` (filter) are **not the same quantity**, and swapping
one for the other silently changes results. The caller's is
`min(|pos − PAM|, |pos + len(ref) − 1 − PAM|)`; the filter's is derived from `indel_info`, measured
from the anchor base only, and therefore always larger (`find_edited_reads.py:2263` documents this).
Measured on the 32-sample cohort they disagree on **162 of 1,498** gated rows — rule 2 drops 838 vs
682 — while both retain all 64 confirmed edits.

The caller's metric is the more correct one, but the 10 bp threshold was calibrated against the
other. So it is exposed as `--cut-dist-source {indel_info,caller}`, defaulting to the calibrated
one, pending recalibration.

### What deliberately did NOT move into the caller

- **The beta-binomial.** It reads `indel_reads`, `total_reads`, `control_indel_reads`,
  `control_reads` — all already columns. It never opens a BAM, so putting it in the caller buys no
  avoided work, forces a global-prior fit into a streaming writer (`find_edited_reads.py` writes
  rows as it goes, `flush=True`), and freezes `bg_rate`/`AQ` behind a 1.5 h rerun. Verified
  separately that the prior is not cohort-dependent: refitting per sample gives **0 of 1,498**
  decision flips at AQ≥5 and an identical queue.
- **The repeat and panel intersects.** `RepeatIndex` loads 5.28M intervals and the panel is a 1 GB
  stream, both currently paid **once per cohort**; per sample they become 32×. `bin/subset_noise_panel.py`
  makes the panel side cheap enough to move later (once per *guide*, not per sample) if wanted.

## Reproduce

```
# Part A, arm 4 (rule 6 standing alone)
bin/review_filter.py results_cart_ponfix/*/*.offtarget_analysis.tsv \
    --repeats rmsk.no_simple.bed GRCh38_no_alt.trf.bed \
    --max-guides 9999 --snv-noise IDPF_...snv.bed.gz \
    --min-reads 2 --min-vaf 0.005 -o arm4.tsv

# Part B
bin/noise_model.py results_cart_ponfix/*/*.offtarget_analysis.tsv --baseline matched -o scored.tsv \
    --truth-wgs "Manual Indel Review/cart_wgs/cart_wgs_merged.xlsx" \
    --truth-ecs "Manual Indel Review/cart_ecs/cart_ecs_merged.csv.gz"
```

## Not done

- `noise_model.py` is **offline**; the `np.select` chain in `review_filter.py` is unchanged.
- No sub-2.5% truth exists in either review, so behaviour below that VAF is still unmeasured.
- A custom baseline built with DRAGEN natively (`--build-sys-noise-vcfs-list`,
  `--build-sys-noise-method=max`) would give the panel format from our own controls; not attempted.
