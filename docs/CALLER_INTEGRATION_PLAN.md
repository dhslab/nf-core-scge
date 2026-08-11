# Caller integration and breakend filtering — plan

Four pieces of work, arising from the review of the 32-sample CAR-T WGS cohort. They are listed in
dependency order; item 2 is the one that unblocks the most.

| # | Goal | Effort | Blocked by |
|---|---|---|---|
| 1 | Move review rules 1–3 into `find_edited_reads.py` | medium | — |
| 2 | Add a `control_bnd_reads` column | small | — |
| 3 | Breakpoint promiscuity rule for BNDs | small | 2 (to be useful) |
| 4 | Use the DRAGEN systematic-noise panels | medium | — |

Everything below was measured on `results_cart_full/` (32 samples, median depth 205×) unless
stated otherwise.

---

## 1. Move rules 1–3 into `find_edited_reads.py`

### Why

Two of the three already exist there in some form, and the duplication is a liability — the caller
and the filter can disagree about the same concept.

| Review rule | Already in the caller? |
|---|---|
| 1. matched control clean | Partly. `-x/--max-in-control` (default 0) exists, but it *filters* rather than reports |
| 2. within N bp of a PAM | **Yes** — `-d/--max-mutation-distance`, default **25**. `review_filter.py` re-applies it at 10 |
| 3. ≥3 distinct indel lengths | No |

Rule 2 is pure duplication: the caller already discards anything beyond 25 bp, then the filter
discards again at 10 bp. One knob, applied twice, in two files.

### The `-x` problem, and why it matters here

`-x/--max-in-control` defaults to `0` and **removes** any event with control support before the
event list is written. The consequence is that the per-event `control_alt_counts` field (field 11
of `indel_info` and `bnd_info`) is **0 on every event, always** — verified across all 2,802 indel
events and 892 BND events in the cohort. It carries no information.

The meaningful statistic is the site-level `control_indel_reads` column, which is computed
independently and is nonzero on 313 rows. That is the column the reporting fix repaired, and it is
what rule 1 actually reads.

**Goal:** make the caller *report* control support rather than silently filter on it. Keep `-x` as
an opt-in, default it to off, and let downstream decide. A filter that erases its own evidence
cannot be audited.

### Depth sensitivity — do not move rule 3 naively

Rule 3 requires ≥3 distinct indel lengths at a site. **A site with 2 indel reads can produce at
most 2 distinct lengths, so it can never pass.** Measured: of the 837 cohort rows with exactly 2
indel reads, the maximum `n_distinct_len` observed is 2, and **zero** pass rule 3.

This makes rule 3 a *de facto* read-depth threshold wearing a different hat. It is safe at
`indel_reads ≥ 10` and meaningless below about 5. If it moves into the caller it must either carry
an explicit minimum-read guard or be expressed depth-independently — for example as a normalised
length-entropy rather than a raw distinct count.

### What should *not* move

Rules 4 and 5 must stay downstream. Both need the whole cohort:

- The panel of normals is built by a separate pass over **all** unedited samples. `find_edited_reads.py`
  runs per-sample, so hosting the PoN there means a two-pass design inside a single-pass script.
- The repeat mask is a static annotation and is cheap to apply once over a merged table.

### Acceptance

- `-d` is read from one place and applied once.
- `control_indel_reads` and `control_bnd_reads` are always populated regardless of `-x`.
- Rule 3 either carries a read-count guard or is replaced by a depth-normalised statistic, with the
  cohort re-scored to confirm the 81-site queue is unchanged.

---

## 2. Add a `control_bnd_reads` column

### Why

There is no site-level control statistic for breakends. `control_indel_reads` exists;
**`control_bnd_reads` does not.** Rule 1 — the workhorse, responsible for 49 of 79 drops (62%) on
indels — therefore has **no breakend equivalent at all**.

Every one of the 2,403 BND events in the cohort reports control support of 0, which is an artifact
of `-x`, not a measurement. We currently cannot distinguish a breakend unique to the edited sample
from one present in the matched control.

### Approach

Mirror the existing indel path: count reads in the matched control that support the breakend
junction, at the site level, computed **before** any `-x` filtering. The junction has two ends, so
the natural statistic is reads spanning or clipped at either breakpoint that also carry the mate
signature.

### Acceptance

- The column is present and nonzero on at least some sites in a cohort known to contain germline SVs.
- Applied to the 32-sample cohort, it reduces the 2,403 BND events measurably. Any reduction is
  informative; today the rule cannot fire at all.

**This is the highest-value item on the list.** It is a small change that turns the single most
effective rule on for an entire event class.

---

## 3. Breakpoint promiscuity rule

### Why

Breakend support is thin — median **1** supporting read across the cohort, only 36 of 2,403 events
have ≥3 — so read count alone does not separate signal from noise. But artifacts have a distinctive
shape: **one breakpoint partnering with many unrelated loci**, which is the signature of an
alignment hub rather than a real junction.

Measured: binning breakpoints at 1 kb and counting distinct partner loci cohort-wide,
`chr1:246,009,982` has **45 distinct partners**. Ten bins of 1,718 have ≥5.

### The rule

```
n_partners(bin) = |{ distinct partner loci for this 1 kb breakpoint bin }|
promiscuous     = n_partners >= 5
```

This is the breakend analogue of rule 3: not "do the reads agree on a length" but "does this
breakpoint agree on a partner".

### Measured effect

```
2,403 BND events
  → reads ≥ 3                    36
  → within 10 bp of a PAM        29
  → not promiscuous (<5)         25
```

It drops exactly the 4 interchromosomal CTLA4 artifacts — all of which share the
`chr1:246,009,982` hub — and nothing else.

### What the 25 survivors are

All same-chromosome, all on-target-to-on-target, in 8 samples (ARID4A, BRAF, IKZF2, KLF12, PDCD4,
PDE7A, RXRB, ZEB2), spanning **413 bp to 126 kb**, median 3 kb.

These are **multi-cut deletions**: two cuts from the same guide's target set, with the intervening
segment excised. The pipeline has been emitting them all along without interpreting them. They
should be reported as an editing outcome in their own right, not as off-target candidates.

### Acceptance

- Promiscuity is computed per run and written alongside the BND table.
- On-target-to-on-target junctions are classified as multi-cut deletions and reported separately
  from off-target breakend candidates, the same way on-target sites are exempted from rules 4 and 5.

---

## 4. DRAGEN systematic-noise panels

Location: `/storage2/fs1/dspencer/Active/shared/refdata/hg38/dragenfiles`

### SV noise — measured coverage of our breakends

| File | Records | Either-end hit @200 bp slop |
|---|---|---|
| `WGS_hg38_v3.1.0_systematic_noise.sv.bedpe.gz` | 311,395 | 13.4% |
| `IDPF_WGS_hg38_v3.0.0_systematic_noise.sv.bedpe.gz` | 2,626,364 | **89.6%** |
| `WGS_FF_Heme_hg38_v3.1.0_systematic_noise.sv.bedpe.gz` | 2,195,842 | 80.4% |

BEDPE columns: `chrom1 start1 end1 chrom2 start2 end2 name score strand1 strand2 precise|imprecise`

### The trap

**Do not apply IDPF naively.** At 200 bp slop it flags 90% of all breakends — including **all 25 of
the real multi-cut deletions**. Used as a blacklist it would delete the entire finding.

The small `WGS_hg38_v3.1.0` panel is the specific one. Any use of these files needs:

- tight slop (test 0 and 50 bp before 200),
- both-ends matching rather than either-end,
- an on-target exemption, exactly as rules 4 and 5 have.

### Goal

Treat the DRAGEN panels as a **supplement** to the run's own PoN, never a replacement. The
self-built PoN is guide-matched; these are generic. The existing `pon_coverage` guard already
distinguishes "this panel says clean" from "this panel does not cover these coordinates", and the
same discipline applies here.

### Open question — the SNV noise BED

`IDPF_WGS_hg38_v.2.0.0_systematic_noise.snv.bed.gz` (1 GB) is untested. It is an SNV panel, so the
prior is that it will behave like the GATK 1000G PoN did — flagging real edits alongside artifacts,
because an SNV-artifact map is not an indel-artifact map. That test was already run once for
GATK and it failed on exactly those grounds.

**Investigate, but expect a negative result**, and measure it against real edits before adopting:
the acceptance bar is that it flags artifacts *without* flagging any of the 77 known on-target
sites.

---

## Summary of what is worth doing first

1. **`control_bnd_reads`** — small change, turns on the most effective rule for a whole event class.
2. **Promiscuity + multi-cut deletion reporting** — no new data needed, and it surfaces a finding
   already sitting in the output.
3. **Caller consolidation of rules 1–3** — real cleanup, but carries the rule-3 depth trap.
4. **DRAGEN panels** — useful as a supplement, dangerous as a blacklist, and the SNV BED is likely
   a dead end worth ruling out cheaply.
