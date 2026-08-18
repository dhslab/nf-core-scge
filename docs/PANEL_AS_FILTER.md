# The DRAGEN systematic-noise panels as filters

Two asks, answered by measurement:

1. Use `IDPF_WGS_hg38_v.2.0.0_systematic_noise.snv.bed.gz` as a filter and count how many indel
   calls overlap a noisy SNV locus.
2. The same for breakends, with `IDPF_WGS_hg38_v3.0.0_systematic_noise.sv.bedpe.gz`.

And behind both, the standing question: **can the panel replace the matched control?**

Everything here is reproducible with [`bin/panel_overlap.py`](../bin/panel_overlap.py), which also
restores a lost capability — the original three-panel measurement in
[`CALLER_INTEGRATION_LOG.md`](CALLER_INTEGRATION_LOG.md) had **no surviving generator**, so nothing
could be re-measured when the caller changed underneath it.

```bash
python3 bin/panel_overlap.py \
  --tables 'results_cart_bnd/*/*.offtarget_analysis.tsv' \
  --queue-all results_cart_bnd/review/review_queue_all.tsv \
  --snv-noise .../IDPF_WGS_hg38_v.2.0.0_systematic_noise.snv.bed.gz \
  --truth-wgs '.../cart_wgs_merged.xlsx' \
  --sv-noise .../WGS_hg38_v3.1.0_systematic_noise.sv.bedpe.gz \
  --sv-noise .../IDPF_WGS_hg38_v3.0.0_systematic_noise.sv.bedpe.gz \
  --sv-noise .../WGS_FF_Heme_hg38_v3.1.0_systematic_noise.sv.bedpe.gz
```

---

## The short answer

**No, and the reason is not panel quality — it is what kind of prior each instrument can supply.**

- **Germline needs a *sample-specific* prior.** A donor's genotype is a property of that donor.
  The correct prior for "is there a real allele here in *this* person" can only come from *this*
  person's own unedited material. A 46-donor panel reports a population average and is
  structurally incapable of representing one donor's genotype, at any panel size.
- **Machine error needs a *population* prior — and here the panel is the better instrument.**
  Sub-detection error rate is donor-independent physics. One control at 157× cannot resolve it;
  46 donors pooled can. See [§7 of the validation doc](NOISE_MODEL_VALIDATION.md#7-what-to-actually-do).

Both statements are true at once, and that is the whole answer.

---

## Ask #1 — the SNV panel as an indel filter

Panel: `IDPF_WGS_hg38_v.2.0.0_systematic_noise.snv.bed.gz` (~1 GB, 46 donors per its own
`##PON SAMPLES:` header). Matching at ±2 bp, as the pipeline does.

Despite "snv" in the filename this panel is **not** SNV-only — its allele column carries `D`
(1,299,670 records) and `I` (892,925) codes. Three nested definitions of "overlap" are therefore
worth separating:

| population | n | any panel record | indel-capable record | shipped rule (≥3 donors & D/I) |
|---|---:|---:|---:|---:|
| all site-rows | 99,238 | 11,757 (11.9%) | 152 (0.15%) | 19 (0.02%) |
| gated rows | 479 | 126 (26.3%) | 46 (9.6%) | 8 (1.7%) |
| **review queue** | **86** | **1 (1.2%)** | **0** | **0** |

**The headline number the ask wanted: 152 of 99,238 indel calls (0.15%) sit on an indel-capable
noisy locus, and only 19 clear the shipped ≥3-donor bar.** Restricting to the rows a human actually
reviews, the queue is essentially panel-clean — one row touches any panel record at all, and none
touches an indel record.

> **Loader check.** The shipped `load_snv_noise` keeps only ONE record per position (the one with
> the most donors), so an SNV-only record with 20 donors could shadow an indel record with 2 and
> make the rule silently under-flag. Measured directly: **19 sites either way — no record is
> shadowed at any queried position on this run.** The lossy retention is not costing anything here,
> but it is a latent defect worth knowing about.

### Against the curated truth label

Scored as a *standalone* filter — flagging means "call this an artifact".

> **Two populations, and it matters which one is quoted.** The label can be joined against every
> site-row the caller emitted, or only against the rows that clear the gate. They answer different
> questions and the panel looks materially different on each, so both are given. **`review_filter.py`
> only ever runs this rule on gated rows**, so the second table is the operational one; the first
> answers the broader "could this replace the noise model outright?". Note that
> [`NOISE_MODEL.md` §6](NOISE_MODEL.md#what-it-is-worth-in-practice--reported-honestly) uses the
> gated population (135 rows) for the cut-distance work — quoting a sensitivity from one table
> against the other is the easy mistake here.

**All site-rows** — 279 rows join the label (69 confirmed / 210 human-rejected):

| definition | flagged | sensitivity | specificity | precision | real edits flagged |
|---|---:|---:|---:|---:|---:|
| *(flag everything)* | 279 | 1.000 | 0.000 | 0.753 | 69 |
| any panel record | 105 | 0.495 | 0.986 | 0.990 | **1** |
| indel-capable record | 36 | 0.171 | 1.000 | 1.000 | **0** |
| shipped rule (≥3 donors & D/I) | 7 | 0.033 | 1.000 | 1.000 | **0** |

**Gated rows only — the operational population** — 135 rows join (64 confirmed / 71 rejected):

| definition | flagged | sensitivity | specificity | precision | real edits flagged |
|---|---:|---:|---:|---:|---:|
| *(flag everything)* | 135 | 1.000 | 0.000 | 0.526 | 64 |
| any panel record | 53 | **0.732** | 0.984 | 0.981 | **1** |
| indel-capable record | 31 | **0.437** | 1.000 | 1.000 | **0** |
| shipped rule (≥3 donors & D/I) | 3 | 0.042 | 1.000 | 1.000 | **0** |

**The panel is ~1.5× more sensitive on gated rows** (0.73 vs 0.50 at the loosest definition, 0.44
vs 0.17 indel-capable). That is not the panel improving; it is the denominator changing. Gated rows
are already enriched for the recurrent-locus artifacts a population panel can see, whereas the full
table is dominated by sites carrying no indel evidence at all — artifacts the panel was never going
to catch. Read the base-rate row first in both cases: at 75% negative the all-rows table hands a
do-nothing filter precision 0.753, and at 53% negative the gated table hands it 0.526.

The conclusions are the same under either denominator, which is the point:

- **The panel is very safe and not very powerful.** Even the loosest definition flags exactly **1**
  confirmed edit in both populations; the two indel-aware definitions flag **none** in both. But
  sensitivity tops out at 0.73 even where the panel is strongest.
- **The rule that actually ships catches almost nothing** — 0.033 and 0.042. On this cohort it
  claims **zero** queue rows.
- **It is a specificity instrument, not a sensitivity one.** That is a perfectly good thing to be —
  it is why the rule sits last in `review_filter.py`, where it can only claim rows the other five
  rules kept — but it means the panel cannot carry the filtering load on its own.

For context, the previously published baseline stands: only **54 of 1,498 gated loci (3.6%)** had
an indel-capable record, and a panel-only noise model scored precision 0.312 against 0.681 for the
matched control ([`NOISE_MODEL_EXPERIMENT.md:62-76`](NOISE_MODEL_EXPERIMENT.md)).

### Can the panel stand in for rule 1?

Rule 1 — "germline, present in the matched control" — is the single biggest filter in the stack,
dropping **154 of 479** gated rows. If the panel could reproduce those drops, the matched control
would be replaceable. It cannot:

```
rows rule 1 dropped as germline : 154
  panel flags, any record       :  55   (36%)
  panel flags, indel-capable    :  34   (22%)
  panel flags, shipped rule     :   6   ( 4%)

germline rows invisible to the panel entirely : 99/154 = 64%
donor support where the panel DOES see them   : median 1, max 35 (of 46)
```

**64% of the germline calls are invisible to the panel**, and of the 36% it does touch, the typical
record is supported by a **single donor** — i.e. a coincidence, not a population statement. The
shipped rule recovers 4%.

This is exactly what the category distinction predicts. A population panel can only recognise a
population-*common* allele; a donor's private variants are, by definition, not in it. Growing the
panel does not fix this, because the thing being asked for is not in the panel's sampling frame.

---

## Ask #2 — the SV panels as a breakend filter

Three panels, on the 1,022 junctions the caller emitted and the 25 that clear `reads ≥ 3`.

> **The on-target exemption is disabled for this measurement.** The shipped rule exempts
> `is_target == 1` junctions from the panel check, and **all 25 gated junctions on this cohort are
> `is_target == 1`** — so with the exemption on, the answer would be trivially zero and would say
> nothing about the panels. The question asked was about panel discrimination, so the exemption is
> turned off and the numbers below are for the panels themselves.

| panel | BEDPE records | merged intervals | **genome covered** |
|---|---:|---:|---:|
| `WGS_hg38_v3.1.0` | 311,395 | 302,943 | **33.7 Mb — 1.09%** |
| `IDPF_WGS_v3.0.0` | 2,626,364 | 998,548 | **1,942 Mb — 62.66%** |
| `WGS_FF_Heme_v3.1.0` | 2,195,842 | 1,375,249 | **1,554 Mb — 50.12%** |

Flagged counts (either end in a panel interval):

| panel | slop | all 1,022 junctions | 25 gated | artifacts (n=2) |
|---|---:|---:|---:|---:|
| `WGS_hg38_v3.1.0` | 0 | 58 (5.7%) | **0 / 25** | 2 / 2 |
| | 50 | 77 (7.5%) | **0 / 25** | 2 / 2 |
| | 200 | 120 (11.7%) | 3 / 25 | 2 / 2 |
| `IDPF_WGS_v3.0.0` | 0 | 852 (83.4%) | 25 / 25 | 2 / 2 |
| | 50 | 878 (85.9%) | 25 / 25 | 2 / 2 |
| | 200 | 934 (91.4%) | 25 / 25 | 2 / 2 |
| `WGS_FF_Heme_v3.1.0` | 0 | 709 (69.4%) | 22 / 25 | 2 / 2 |
| | 50 | 745 (72.9%) | 25 / 25 | 2 / 2 |
| | 200 | 835 (81.7%) | 25 / 25 | 2 / 2 |

### The mechanism, which is simpler than previously stated

**IDPF v3.0.0 covers 62.66% of hg38.** A hit against it is close to uninformative by construction —
it flags 83% of *all* junctions, real and artifactual alike. `WGS_FF_Heme` covers 50% and behaves
the same way. Only `WGS_hg38_v3.1.0`, at **1.09% coverage**, is selective enough for a hit to mean
anything, and it is the only panel that discriminates: **0 of 25 real junctions, 2 of 2 available
artifacts**, at slop 0 and 50.

This is a sharper explanation than the one in `CALLER_INTEGRATION_LOG.md`, which framed IDPF as
"actively harmful" because it flags real junctions at a higher rate than noise. That framing is
still arithmetically true, but the cause is not a population mismatch — it is that the panel's
intervals blanket most of the genome, so almost everything hits.

### Would a `min_donors` rule rescue the large panels? No.

The BEDPE carries **no donor-count column and no `##PON SAMPLES:` header**, so there is no direct
analogue of the SNV panel's `n_donors` — the very statistic that makes the SNV panel usable. But it
is recoverable: field 7 is a candidate name that embeds the donor it came from
(`ImpreciseNoiseCandidate_LP7108672-DNA_A06_42_DRAGEN:BND:...`), so pooling records by interval and
counting distinct donors reconstructs the missing column.

Donor support behind the hits on the 25 gated junctions (slop 50):

| panel | hit | 1 donor | ≥2 donors | ≥3 donors | max |
|---|---:|---:|---:|---:|---:|
| `WGS_hg38_v3.1.0` | 0 | 0 | 0 | 0 | 0 |
| `IDPF_WGS_v3.0.0` | 25 | **0** | 25 | **25** | 17 |
| `WGS_FF_Heme_v3.1.0` | 25 | 6 | 19 | 12 | 9 |

**No threshold rescues IDPF.** Every one of its hits on a real junction is supported by ≥3 donors
(up to 17), so these are genuinely recurrent loci in its baseline, not one-donor coincidences. A
`min_donors ≥ 3` rule would leave IDPF at 25/25 and take FF_Heme only from 25/25 to 12/25 — still
unusable. The panel *choice* is load-bearing on the SV side in a way it is not on the SNV side, and
no tuning substitutes for picking the right file.

Also worth noting: matching is a bare either-end interval hit
([`review_filter_bnd.py:130-135`](../bin/review_filter_bnd.py)) — there is no requirement that the
*partner* end match too, which would be a much stronger test and is not currently available.

### ⚠️ Flagged for a decision, not changed here

`params.sv_noisefile` ([`nextflow.config:208`](../nextflow.config)) is **IDPF v3.0.0**, and it is
passed to DRAGEN itself as `--sv-systematic-noise`. That is the panel measured above as covering
63% of the genome and flagging 100% of the real junctions on the review side.

The review-side default is already the safe panel (`review_sv_noise` in `logs/run_cart_bnd.sh`
points at `WGS_hg38_v3.1.0`), so this is about what DRAGEN does upstream, which is a different
question and out of scope here. **It is worth a decision.**

### Reproducibility note — an honest discrepancy

`CALLER_INTEGRATION_LOG.md:587-595` published this table against **36 gated junctions (25 real, 11
artifacts)** from the 2,403-junction set of 2026-08-10:

| panel | published: real flagged | published: artifacts | **this run: real flagged** |
|---|---|---|---|
| `WGS_hg38_v3.1.0` | 0 / 25 | 11 / 11 | **0 / 25** ✓ |
| `IDPF_WGS_v3.0.0` | 25 / 25 | 11 / 11 | **25 / 25** ✓ |
| `WGS_FF_Heme_v3.1.0` | 24 / 25 | 11 / 11 | **25 / 25** ✗ |

- Raw record counts reproduce exactly (311,395 / 2,626,364 / 2,195,842).
- Two of three "real flagged" figures reproduce exactly.
- **FF_Heme differs by one** (24 → 25). The junction set changed (2,403 → 1,022) and the caller
  changed with it, so this is not the same 25 junctions. Reported, not reconciled.
- **The 11-artifact set no longer exists on disk.** No run retained it: the current run drops
  *zero* junctions to rules 2–4, so it has no artifacts at all. The only recoverable artifact set
  is 2 junctions from `results_cart_noise`, and all three panels flag both. **The artifact column
  of the published table is therefore not currently reproducible**, and any future claim resting on
  "11/11" should be re-derived rather than quoted.
