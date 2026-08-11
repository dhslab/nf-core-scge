# Lowering the evidence gate — expected impact

**Change:** `indel_reads >= 10 & indel_fraction >= 0.05` → `indel_reads >= 2 & indel_fraction >= 0.005`
(`review_min_reads`, `review_min_vaf` in `nextflow.config`).

Every number below is **measured**, not projected: `review_filter.py` was run directly over the same
32 `*.offtarget_analysis.tsv` inputs, the same panel of normals and the same repeat BEDs, with only
the two thresholds changed.

> **Confirmed by the pipeline.** SLURM job 2707025 (`-profile ris2,apptainer`, outdir
> `results_cart_lowgate/`) completed successfully, 2026-08-10. It produced **1,498 gated rows and a
> 91-site queue**, and the queue is **set-identical** to the simulation below — same samples, same
> coordinates, no losses against the previous 81. 91 snapshots rendered. The old run in
> `results_cart_full/` is untouched.

---

## Headline

| | old gate | new gate | |
|---|---|---|---|
| Rows clearing the gate | 160 | **1,498** | 9.4× |
| Review queue | 81 | **91** | +10 |
| On-target | 77 | 81 | +4 |
| **Off-target candidates** | **4** | **10** | **2.5×** |
| Sites lost | — | **0** | all 81 retained |

**The queue grows 12%; the work grows 150%.** The 81 original sites all survive — no regression —
but the number of rows that actually need a human decision goes from 4 to 10, because on-target
sites are confirmations rather than decisions. That is the number to watch.

### Per-rule drops

| Rule | old | new |
|---|---|---|
| germline (in control) | 49 | 352 |
| far from PAM | 18 | **680** |
| single indel length | 9 | **359** |
| known-bad (PoN) | 1 | 6 |
| repeat region | 2 | 10 |
| **kept** | **81** | **91** |

The five rules absorb 1,407 of the 1,498 gated rows. Rule 2 becomes the workhorse — at low VAF most
new material is scattered background indels that are simply nowhere near a cut site.

---

## Finding 1 — the `indel_reads >= 2` half of the change is inert

**Rule 3 requires ≥3 distinct indel lengths. A site with 2 indel reads carries at most 2 distinct
lengths, so it can never pass.**

Measured: **837 rows** enter with exactly 2 indel reads. Maximum `n_distinct_len` observed among
them is 2. **Zero** survive.

| indel_reads | rows | can satisfy rule 3? |
|---|---|---|
| 2 | 837 | **no — arithmetically impossible** |
| 3 | 225 | only if all 3 reads differ |
| 4+ | 272 | yes |

So 56% of the newly admitted rows are dropped by arithmetic rather than by evidence. Every one of
the 10 new queue sites has **3–10** indel reads: they came from relaxing the *VAF* term, not the
read term.

**Implication:** `indel_reads >= 3` would give an identical queue with 837 fewer rows to carry,
audit and store. If the intent of `>= 2` is genuinely to see 2-read sites, then rule 3 has to change
at the same time — it is currently acting as a hidden `indel_reads >= 3` filter.

---

## Finding 2 — a new asymmetry between the gate and the panel of normals

Lowering the review gate without lowering the PoN threshold makes the panel **systematically
blinder than the filter it feeds**.

The PoN blacklists a site needing `>=3 reads AND >=2% VAF` in an unedited sample (`pon_min_reads`,
`pon_min_vaf`). The review gate now admits sites at 2 reads. So there is a band — 2 reads in the
control — where the filter can see a site but the panel structurally cannot flag it.

**This is not hypothetical; it already leaks.** `chr18:61,959,927` appears in **28 of 32 samples
under 23 distinct guides** — the textbook definition of a recurrent artifact. Its PoN record:

```
chrom   pos        n_donors  max_indel_reads  max_indel_fraction  blacklisted
chr18   61959928   0         2                0.0526              0
```

The controls **do** carry it — 2 reads, 5.3% VAF — but `pon_min_reads = 3` misses it by a single
read, so it is not blacklisted. Under the old gate it never cleared and the gap was invisible.
Under the new gate **3 copies reach the queue** (CART_NS0027-B2M_1, IKZF2, PLCB2), and they are
3 of the 10 off-target candidates.

### Two independent fixes, either of which closes it

1. **Lower `pon_min_reads` 3 → 2** to match the gate. chr18 then has `max_indel_reads = 2 ≥ 2` and
   `max_indel_fraction = 0.0526 ≥ 0.02`, so it blacklists. **Keeping the two thresholds in step is
   the general principle**: the panel must be at least as sensitive as the filter consuming it.
2. **Apply cross-guide recurrence *alongside* the PoN, not instead of it.** `n_guides = 23` is
   already computed and sitting in the output — it is simply unused whenever a PoN is supplied
   (`review_filter.py:224` vs `:228`, an if/else). Making rule 4
   `is_off & (in_blacklist | n_guides >= max_guides)` drops all 3 chr18 rows.

Either takes the queue **91 → 88** and off-target candidates **10 → 7**.

---

## Finding 3 — what the 10 new sites actually are

| Sample | Site | Reads | VAF | On-target | Verdict |
|---|---|---|---|---|---|
| CART_NS0027-CTLA4_2 | chr7:142,792,021 | 8 | 0.096 | ✓ | TRBC1, low-efficiency replicate |
| PLCB2 | chr18:61,959,927 | 8 | 0.079 | | **recurrent artifact** (23 guides) |
| IKZF2 | chr18:61,959,927 | 5 | 0.046 | | **recurrent artifact** |
| CART_NS0027-PD1_1 | chr2:241,852,750 | 8 | 0.040 | ✓ | |
| ZEB2 | chr4:26,561,342 | 6 | 0.040 | | new candidate |
| CART_NS0027-CTLA4_2 | chr2:203,870,831 | 10 | 0.039 | ✓ | CTLA4 |
| CART_NS0027-B2M_2 | chr15:44,711,599 | 5 | 0.023 | ✓ | **the second B2M cut site** |
| CART_NS0027-B2M_1 | chr18:61,959,927 | 3 | 0.021 | | **recurrent artifact** |
| CART_NS0027-CTLA4_1 | chr2:160,064,217 | 4 | 0.018 | | new candidate |
| PLCB2 | chr5:159,263,654 | 3 | 0.015 | | new candidate |

**Two are genuinely valuable.** `chr15:44,711,599` in CART_NS0027-B2M_2 is the *second* B2M cut
site — the deck's TRAC/TRBC slide shows editing split across two candidate B2M positions, and this
change recovers the low-efficiency one, confirming both fired in that sample. The CTLA4 and PD1
on-target additions similarly fill in low-efficiency replicates.

**Three are the chr18 artifact** described above.

**Three are new off-target candidates** (ZEB2 chr4, CTLA4 chr2, PLCB2 chr5) at 1.5–4.0% VAF, each
seen in a single sample. At 205× these rest on 3–6 reads. They are exactly the population the
system cannot currently adjudicate — real enough to demand a look, thin enough that a look will
not settle them.

---

## Operational cost

| | old | new |
|---|---|---|
| `review_queue_all.tsv` rows | 160 | 1,498 |
| Snapshots rendered | 81 | 91 |
| Runtime | ~9.5 min | ~unchanged |

Snapshot rendering scales with the *queue*, not the gated set, so cost barely moves. The audit file
grows ~9×, which is still trivially small.

---

## Recommendation

1. **Use `indel_reads >= 3`, not 2.** Identical queue, 837 fewer inert rows. Going to 2 is only
   meaningful if rule 3 is reworked at the same time.
2. **Lower `pon_min_reads` to match whatever the gate becomes.** The panel must never be less
   sensitive than the filter reading it.
3. **Make rule 4 the union of PoN and cross-guide recurrence.** The signal is already computed and
   currently thrown away, and it catches exactly the class the PoN's read threshold misses.
4. **Do not present the 10 off-target candidates as findings** until 1–3 are applied. Three of the
   ten are one artifact counted three times, which would be a bad slide.

With all three applied the queue is **88** — 81 on-target, 7 off-target — and every one of the
7 is a single-sample, non-recurrent, control-clean candidate.
