# Caller integration — working log

Running notes for the work described in [`CALLER_INTEGRATION_PLAN.md`](CALLER_INTEGRATION_PLAN.md).
Newest entries at the bottom. Every claim here should be reproducible from a command in the entry.

## Safety net

| | |
|---|---|
| Restore point | commit `e9b81b9`, `bin/find_edited_reads.py` **tracked and clean** |
| Undo a bad edit | `git checkout bin/find_edited_reads.py` |
| Uncommitted before we start | `nextflow.config` (gate 2/0.005, `pon_min_reads` 2), `docs/CALLER_INTEGRATION_PLAN.md`, `docs/THRESHOLD_CHANGE_IMPACT.md` |
| Regression baselines | `results_cart_full/` (gate 10/0.05), `results_cart_lowgate/` (gate 2/0.005), `results_cart_ponfix/` (+ `pon_min_reads` 2) |

**Rule for this work:** no edit to `find_edited_reads.py` lands without (a) a baseline capture before,
(b) a diff of the *outputs* after, not just an exit code. The caller is on the shared analysis path —
a regression here breaks the whole pipeline, not just the off-target arm.

## Order of work

Following the plan's own priority:

1. **Item 2 — `control_bnd_reads`** (in progress) — smallest surface, unblocks item 3
2. Item 3 — breakpoint promiscuity
3. Item 1 — consolidate rules 1–3 into the caller
4. Item 4 — DRAGEN noise panels

---

## 2026-08-10 — Item 2 reconnaissance

Goal: add a site-level `control_bnd_reads` column so rule 1 (matched control clean) can apply to
breakends. Today there is `control_indel_reads` but no BND equivalent, so rule 1 — which removes 352
of 1,410 dropped rows on indels — cannot fire on breakends at all.

### What we already know going in

- Per-event `control_alt_counts` (field 11 of both `indel_info` and `bnd_info`) is **0 on every
  event, always**, because `-x/--max-in-control` defaults to 0 and *removes* events with control
  support before they are written. Verified across 2,802 indel and 892 BND events.
- The meaningful statistic is the site-level `control_indel_reads` column, nonzero on 313 rows.
- 2,403 BND events in the 32-sample cohort; median 1 supporting read; 51% interchromosomal.

### Safety net established

**Test suite is runnable and green: `30 passed`.** It was not runnable as documented — `pytest` is
absent from every container on this system (`docker-scge-offtarget`, `docker-scge`,
`docker-baseimage`) and from the system Python. Recipe that works:

```bash
python3 -m venv --system-site-packages pyt      # inherits pandas/pysam/sklearn/scipy/joblib
./pyt/bin/pip install pytest edlib pyranges     # edlib + pyranges are the only real gaps
./pyt/bin/python -m pytest tests/ -q            # -> 30 passed
```

This matters because `tests/conftest.py` builds a **synthetic aligned-read workspace** and
`tests/test_offtarget_glue.py` actually executes `find_edited_reads.py` against it. So we have a
real regression check on the file we are about to edit, not just an import smoke test.

Baseline recorded **before any edit**: `30 passed, 1 warning`. Any change must keep that at 30.

### Reconnaissance findings

**F1 — the control-counting machinery for breakends already exists and is correct.**
`add_normal_counts` (:828) builds an `altseq` for every event via `generate_contig` (:269), and that
function has explicit BND handling (:291): it parses the bracket notation, fetches the local and
remote flanks, reverse-complements where the orientation demands it, and stitches a **chimeric
junction contig**. Control reads are then counted against that contig by the same three-tier match
used for indels. So control support at a breakend is *already being computed properly* — nothing
needs to be written from scratch.

**F2 — it is computed, then silently folded into `control_indel_reads`.**

```
:2225  control_alt_observed = int(indelcounts['control_alt_counts'].sum())   # ALL event types
:2228  indelcounts = indelcounts[control_alt_counts <= max_in_control ...]   # -x filter
:2238  control_indel_reads = control_alt_observed
:2261  bnds   = indelcounts[(alttype == 'BND') & (ref != '.')]               # split happens
:2262  indels = indelcounts[(alttype != 'BND') & (ref != '.')]               # 36 lines later
```

The sum at :2225 runs over `indelcounts` **before** the indel/BND split, so it spans indels, BNDs
and REF placeholder rows alike. **`control_indel_reads` is therefore mislabelled** — it is "control
support at this site across all event types", not "control indel reads".

The summation position itself is deliberate and correct (comment at :2220): it runs before the `-x`
filter so the reported column survives `-x 0`. That is the fix already made. The event-type split is
the part that was never done.

**F3 — the required change is a split, not new logic.** Sum `control_alt_counts` separately over the
BND rows at :2225 and expose it. Roughly four lines.

**F4 — no downstream consumer reads these tables positionally.** All 18 scripts that read
`*.offtarget_analysis.tsv` use `pd.read_csv(..., sep='\t')` with name-based access. The only
positional reads in `bin/` (`usecols=[0,1,2]`) are against **BED** files in `make_hotspot_vcf.py`
and `targets_csv_to_vcf.py`, which are unaffected. **Appending a column to the header is safe.**

**F5 — evidence the mixing is real, but not yet proof it is numerically material.** 4 cohort sites
have `control_indel_reads > 0` with `indel_count == 0` and `bnd_count > 0`
(e.g. `CART_NS0027-B2M_2 chr9:106,439,687`, control 37/110). That is *consistent* with a BND
contribution, but equally consistent with an indel event that the `-x 0` filter removed after its
support was counted. **Do not claim BND contamination of `control_indel_reads` as measured** — the
structural fact (F2) is certain; the magnitude is not. Splitting the sum is what will measure it.

### Proposed change (not yet applied)

Purely **additive**. `control_indel_reads` keeps its current value bit-for-bit, so no existing
behaviour moves and the review filter is unaffected.

1. `:2225` — alongside `control_alt_observed`, compute `control_bnd_observed` over BND rows only.
2. `:2231-2239` — carry it into a `control_bnd_reads` local.
3. `:1889-1891` — append `control_bnd_reads` to `header_columns`, after `bnd_count`/`bnd_info`.
4. `:2317-2322` and the `fp_log` block at `:2305-2310` — emit it in both row writers, in the same
   position as the header.

Deliberately **not** doing in this step: correcting `control_indel_reads` to be indel-only. That
changes existing numbers and needs its own before/after measurement. Logged as a separate decision.

Risks to watch: the `fp_log` writer and the main writer are separate field lists that must stay in
lockstep with the header, and `--enable-crispr-prediction` appends three more columns after ours.

### Applied — 5 edits

`control_bnd_reads` inserted after `bnd_info` (col index 14 of 17), keeping the BND statistics
together. Edits: header (:1889), the BND-restricted sum (:2225), local init (:2231), assignment
(:2238), and **both** row writers.

### Verification

**V1 — structural.** `py_compile` clean. Header = 17 columns; main writer = 17; `fp_log` writer = 17.
The two writers and the header are independently parsed and counted, so they cannot silently drift.

**V2 — regression suite.** `30 passed`, identical to the pre-edit baseline.

**V3 — A/B on real data, identical inputs.** This is the one that matters. Extracted the HEAD version
(`git show HEAD:bin/find_edited_reads.py`) and ran **both** versions on the same CRAMs and the same
92-site target subset for `KLF12-KO-DNA`:

```
rows: HEAD 92, patched 92, same order
HEAD 16 cols, patched 17 cols, added: ['control_bnd_reads']
ALL 16 PRE-EXISTING COLUMNS BIT-IDENTICAL
```

The change is provably additive on real data.

> **Method note, worth remembering.** The first comparison was patched-output vs the *stored Aug 7
> full-cohort* output, and it showed `bnd_info` differing on 11 rows plus 20 rows failing to join —
> `start` had shifted by 1 bp. That was **entirely an artifact of subsetting the target VCF**:
> adjacent targets merge into different intervals when their neighbours are absent, so the site
> boundaries move. It was not caused by the edit. Never A/B a subset run against a full run —
> re-run the baseline binary on the identical subset.

**V4 — the column reads 0 everywhere, and chasing that found a pre-existing bug.**

`control_bnd_reads` is 0 at **all 92 KLF12 sites** and at **all 204 CTLA4_1 sites** (203 of which
carry breakends, 234 events). That is not because breakends lack control support. It is an
**off-by-one in the position passed to `add_normal_counts`**, which predates this work.

*Evidence, in order.*

1. The control CRAM is not clean at these loci. At `chr1:246,009,900-246,010,060` the **control** has
   174 reads, of which 111 survive the read filters, **37 carry `SA` tags**, and **4 point at chr8** —
   the exact partner of the breakend under test.
2. Called directly with the variant, the counting works: `add_normal_counts` returns
   **`control_alt_counts = 11`**.
3. The variant is not a straw man. The pipeline's own `bnd_info` for that site contains
   `pos=246009982 ref=T alt=T[chr8:94633103[ counts=5 ctrl=0` — byte-identical to the probe.
4. `control_total_counts = 104` from the probe equals the pipeline's `control_reads = 104`, so the
   probe reproduces the pipeline's call exactly.
5. Sweeping the position is decisive:

   | position passed | `control_alt_counts` |
   |---|---|
   | 246,009,982 — 1-based, correct | **11** |
   | 246,009,981 — what the pipeline passes | **0** |
   | 246,009,983 | 3 |

**Root cause.** `indelcounts['pos']` is **0-based** — `bnd_info`/`indel_info` keys are built as
`r['pos']+1` (:2265, :2268) precisely to convert it. But `add_normal_counts` treats its `pos` as
**1-based** (`start_idx = pos - 1`, :839). So every contig is built one base off.

**Why this destroys breakends but only degrades indels.** `generate_contig` puts the junction at an
exact point, so a 1 bp shift means a junction-spanning read is no longer a substring of `altseq` and
tier 2 (`read_seq in v['altseq']`) fails outright. Tier 1 also fails, because it requires
`vcf_dict['pos'] == v['pos']` across the same mismatched convention. Tier 3 (edlib) is the only
survivor, and its accept test is indel-shaped — it branches on `ref_len` vs `alt_len`, where for a
breakend `alt_len` is the length of the *bracket string* (`len("T[chr8:94633103[") == 16`), not a
sequence length. So for BNDs nothing can fire. For indels tier 3 still recovers most cases, which is
why `control_indel_reads` is nonzero on 313 cohort rows despite the same shift.

**Status of this item.** The `control_bnd_reads` column is correctly implemented, verified additive,
and currently reads 0 everywhere because of the above. **The column is not the blocker — the
off-by-one is.**

**Deliberately not fixed in this step.** Correcting the position convention would change
`control_indel_reads` on real data, which feeds rule 1 of the review filter and therefore the queue.
That needs its own before/after measurement on the full cohort, exactly like the threshold change
did. Not measured yet: how much `control_indel_reads` moves once the shift is corrected.

**Next actions, in order.**
1. Quantify the indel-side impact of the off-by-one on a full sample (`control_indel_reads` before
   vs after) before touching the convention.
2. Fix the convention at the call site (:2193) rather than inside `add_normal_counts`, so the
   function keeps its documented 1-based contract.
3. Give tier 3 a breakend-aware accept test, since the `ref_len`/`alt_len` comparison is meaningless
   for bracket notation.
4. Only then re-measure `control_bnd_reads` and decide whether rule 1 is worth applying to BNDs.

---

## 2026-08-10 — Measuring the off-by-one fix. **Verdict: do not ship it yet.**

Built a candidate (`scratchpad/bndverify/cand_fer.py`, one line: `pos = int(row['pos']) + 1` in
`add_normal_counts`) and ran it against the current tree on the identical KLF12 92-site subset and
identical CRAMs.

| column | current | candidate | sites changed |
|---|---|---|---|
| `total_reads` | 15,815 | 16,248 | 19 |
| `indel_reads` | 413 | **846** | 19 |
| `indel_count` | 90 | 123 | 17 |
| `control_indel_reads` | 1,100 | **541** | 40 |
| `control_bnd_reads` | 0 | **2** | 1 |
| rows passing gate (2, 0.005) | 24 | 28 | — |

**The diagnosis was right: the fix does turn on breakend control counting** (0 → nonzero). That
confirms the off-by-one was the blocker, not the new column.

**But the blast radius is far larger than the control columns.** `indel_reads` *doubles*. That is
not a side effect of counting — `total_reads`/`indel_reads`/`indel_count` are summed at :2236-2237,
**after** the `-x` filter at :2228. Change control support and you change which events survive the
filter, so the caller's primary variant calls move. This is not a "control column" fix; it changes
what the caller reports as edited.

**Which version is correct is NOT established.** I tried to settle it on ground truth and could not:

- At `chr4:183,359,487` the current build says 90 control indel reads and the candidate says 0. The
  control genuinely is not clean there — 109 of 147 filtered reads carry a 1 bp deletion, and it
  localises to a single germline event (0-based 183,359,514, anchor key 183,359,514, 88 reads).
- But neither number is obviously right, because the *events called at that site* are a 38 bp
  deletion (`C|DEL38`, key 183,359,514) and a 1 bp deletion at a different anchor (key 183,359,510).
  Neither is the germline deletion the control carries. So "the control has indels here" does not
  decide whether either event is supported.
- The current build's 90 may be spurious matching by a mis-centred contig; the candidate's 0 may be
  correct-but-strict. **Unresolved.**

**Coordinate convention, as far as it is established.** `get_cigar_indel_vcf` sets
`out_dict['pos'] = anchor_ref_pos = current_ref + ref_consumed - 1`, seeded from
`read.reference_start` (pysam, 0-based), and fetches `ref` at that same 0-based coordinate. The
`+1` when building info keys converts to 1-based. `add_normal_counts` and `generate_contig` both do
`start_idx = pos - 1`, i.e. expect 1-based. **The inconsistency is real.** For the breakend case it
was verified end-to-end (correct position → 11 control reads and the control demonstrably carries 4+
reads to the exact partner; the position the pipeline passes → 0). For indels the equivalent
end-to-end verification has *not* been done.

### Decision

**Do not apply the positional fix on the current evidence.** A change that doubles `indel_reads`
must be validated against known truth, not reasoned about. It feeds the evidence gate, the review
queue and every downstream number in the deck.

**Validation set to use next:** the manual-review truth (`Manual Indel Review/`, the
`manual_review` column) and/or the AAVS1 curated truth (`crispr_ml/AAVS1_training_{tp,tn}.tsv`).
Run both builds over a cohort with known calls and compare recall/precision against the labels.
Whichever build agrees with the manual calls is correct. That is a measurable question and should
not be answered by inspection.

### Repository state

`bin/find_edited_reads.py` contains **only** the `control_bnd_reads` addition, which is verified
additive (all 16 pre-existing columns bit-identical). The positional fix exists solely as a
scratchpad candidate and has **not** been applied. Nothing is committed.

---

## 2026-08-10 — The self-consistency test, and a correction to the entry above

The manual-review truth turned out to be unusable for this question: `KLF12-KO-DNA.edited_reads.xlsm`
has 2,328 rows and a `pass` column with **2** annotations, both on-target. It is a positives-only
set — recall only, no negatives, so it cannot adjudicate over- vs under-counting. (Consistent with
what is already known: 55 confirmed positives against tens of thousands of NaNs.)

**A label-free test settles it instead.** Run the caller with `--edited-bam` and `--control-bam` set
to the *same* CRAM. Every event called in "edited" is then, by construction, fully supported in
"control". Any event surviving the `-x 0` filter is a provable miss by the control matcher, and
`control_indel_reads` should be high. No labels needed.

Both builds, KLF12 92-site subset, identical CRAM on both sides:

| metric | current | global `+1` candidate | ideal |
|---|---|---|---|
| events escaping `-x 0` | **12** | 72 | ~0 |
| `indel_reads` | **77** | 465 | ~0 |
| `control_indel_reads` | **3925** | 2000 | high |
| `control_bnd_reads` | 5 | **40** | high |
| sites with ≥1 escape | **9 / 92** | 27 / 92 | 0 |

### Correction to the previous entry

**The earlier conclusion — "`add_normal_counts` receives 0-based positions and treats them as
1-based" — was too broad.** It holds for **breakends**, where it was verified end-to-end. It does
**not** hold for indels: the current handling is substantially *more* correct there, and the global
`+1` makes indels six times worse by this measure. The two producers of `pos` disagree by one base,
and `generate_contig`'s BND branch and linear branch inherit that disagreement:

- BND path (`:395-410`): `bp1_pos = L['r_end']` (0-based), `ref_base` fetched at `bp1_pos`, and the
  code comments `'pos': bp1_pos, # VCF 1-based POS is usually this value`.
- Indel path (`get_cigar_indel_vcf`): `pos = anchor_ref_pos` seeded from `read.reference_start`,
  `ref` fetched at the same coordinate.

They look alike on the surface, which is why reading the code was not enough and the measurement
was. **Do not fix this globally.** The `+1` belongs to breakends only.

### The BND-only fix — confirmed and applied

`if alttype == 'BND': pos += 1`. Self-consistency, all 92 sites:

| metric | current | global `+1` | **BND-only** | ideal |
|---|---|---|---|---|
| `indel_count` (escapes) | 12 | 72 | **12** | ~0 |
| `indel_reads` | 77 | 465 | **54** | ~0 |
| `control_indel_reads` | 3925 | 2000 | **3960** | high |
| `control_bnd_reads` | 5 | 40 | **40** | high |
| `bnd_count` (escapes) | 38 | 22 | **22** | ~0 |
| sites with escaped indel events | 9 | 27 | **9** | 0 |

It takes the best of both: **indel event calling identical to current** (`indel_count` equal
row-by-row, 9 sites with escapes either way) while breakend control support improves **8×** and
breakend escapes drop from 38 to 22.

The cleanest confirmation is the arithmetic: `control_indel_reads` moves by **+35**, and
`control_bnd_reads` gains **+35**. Identical. Since `control_indel_reads` is summed over all event
types (finding **F2** at the top of this log), the *only* thing that changed is the newly-recovered
breakend support. Nothing on the indel side moved at all.

**Applied to `bin/find_edited_reads.py`.** `py_compile` clean, test suite **30 passed**.

### Real-data A/B (edited vs true control), KLF12, 92 sites

| column | before | after | sites changed |
|---|---|---|---|
| `indel_count` | 90 | 90 | **0** |
| `total_reads` | 15,815 | 15,814 | 1 |
| `indel_reads` | 413 | 412 | 1 |
| `control_indel_reads` | 1,100 | 1,102 | 1 |
| `bnd_count` | 41 | 40 | 1 |
| `control_bnd_reads` | 0 | **2** | 1 |

Precisely targeted: **one** breakend at **one** site is now correctly recognised as present in the
matched control and removed by `-x`. Indel event calling is untouched (`indel_count` 90 → 90 with
zero sites changed), and `control_indel_reads` moves by exactly the `control_bnd_reads` gain.

Note the scale honestly: on this sample the change is nearly a no-op, because KLF12's breakends are
single-read junctions at homology sites where the control genuinely carries nothing. The fix matters
where the control *does* carry the junction.

### Where it actually bites — CART_NS0027-CTLA4_1

Re-ran the sample's 204 breakend-bearing sites. At the `chr1:246,009,98x` alignment hub:

| | before | after |
|---|---|---|
| `bnd_count` | **25** | **2** |
| `control_bnd_reads` | **0** | **197** |

**23 of the 25 breakends at that site are now correctly identified as present in the matched control
and removed.** Across all 204 sites:

| column | before | after |
|---|---|---|
| `indel_count` | 42 | **42** (identical row-by-row) |
| `bnd_count` | 233 | **210** |
| `control_bnd_reads` | 0 | **197** |
| `control_indel_reads` | 37 | 234 |

23 spurious breakends removed cohort-wide for this sample, indels untouched. `control_indel_reads`
moves only because it sums over all event types (finding **F2**).

This is the same locus flagged earlier by breakpoint promiscuity (45 distinct partners) and the same
one behind all four interchromosomal "off-target" breakend candidates in the CTLA4 samples. **Two
independent rules — matched control and promiscuity — now reject it.** That is the strongest
possible outcome: the artifact class item 3 was designed to catch is also caught by rule 1 once rule
1 can actually see breakends.

### Why this was not findable by reading the code (item 2)

Both producers look identical on inspection — each stores a 0-based coordinate and fetches `ref` at
it. The disagreement is in `generate_contig`, whose BND and linear branches consume `pos`
differently. Two careful code readings gave the wrong answer here; the self-consistency test gave
the right one in a single run. **For any further coordinate work in this file, measure it — the
edited==control trick costs about three minutes and needs no labels or truth set.**

---

## 2026-08-10 — Item 3: breakpoint promiscuity. Built, and it is non-binding.

New script: **`bin/review_filter_bnd.py`**. Same shape as `review_filter.py` — expands `bnd_info`
into one row per junction, gates, applies rules, writes a queue and an optional `--keep-all` audit.

Cohort result (32 samples, pre-fix inputs):

```
input junctions       : 2404
cleared the gate      :   36   (reads>=3)
  dropped, far from PAM            : 7
  dropped, promiscuous breakpoint  : 4
BREAKEND QUEUE        :   25
    multi-cut deletion    :  23 in 7 sample(s)   span 413-125,747 bp
    deletion at cut site  :   2 in 1 sample(s)
```

**There are no off-target breakend junctions in this cohort.** Everything that survives is an
on-target editing outcome.

### The promiscuity rule does not earn its keep — and that is item 2's doing

Measured, not assumed:

- After rule 1 (`reads >= 3`), the **only** promiscuous breakpoint cohort-wide is `chr1:246,009,98x`.
  Every other breakpoint with ≥5 partners has a maximum of **1** supporting read, so rule 1 removes
  them first.
- That one hub is exactly what item 2's control fix now removes on its own. Re-running three samples
  through the fixed caller collapses the signal:

  | sample | junctions before → after | max partners before → after |
  |---|---|---|
  | CART_NS0027-CTLA4_1 | 234 → 210 | **25 → 3** |
  | CART_NS0027-CTLA4_2 | 237 → 201 | **31 → 2** |
  | CART_NS0027-Regnase1_2 | 192 → 191 | 2 → 2 (no hub to begin with) |

So rule 3 is **kept but currently non-binding**: it costs nothing, it is the only defence against a
hub that happens to be absent from the matched control, and its `n_partners` column is worth
reporting. `0 dropped by rule 3` is the expected reading, not a broken rule. This is documented in
the script's docstring so nobody "fixes" it later.

### The on-target exemption is not optional

**6 of the 10 breakpoints with ≥5 distinct partners are the intended cut sites** — TRAC
(`chr14:22547`, 22 partners), TRBC1/TRBC2 (`chr7:142792`, `chr7:142801`), B2M (`chr15:44711`).
A real Cas9 cut throws junctions everywhere, so the true edits are among the most promiscuous
breakpoints in the cohort. Applying rule 3 without exempting on-target sites deletes the findings.

### Classification, corrected once during development

`is_target` describes the **near** end — the site the junction was found at — so a junction anchored
at a cut site is an on-target outcome even when its far end is unremarkable. The first version
labelled anything without a cut site at *both* ends "off-target junction", which mislabelled the
BRAF pair (`chr7:140,801,460 → chr7:140,834,645`, 33 kb, 7 and 3 reads). BRAF's six on-target sites
do not include the far end, so that is a single cut resected and joined 33 kb out — a real on-target
deletion, not an off-target. Final taxonomy:

| both ends at cut sites, same chrom | multi-cut deletion |
|---|---|
| near end at a cut site, same chrom | deletion at cut site |
| near end at a cut site, other chrom | translocation at cut site |
| near end not a cut site | off-target junction |

> **Superseded 2026-08-20 — these are inversions, not deletions.** The taxonomy above classifies on
> position only and never reads `strands`. Every same-chromosome row in this queue is `+-`/`-+`,
> the inverted adjacency; a deletion is collinear (`++`/`--`) and, on the same chromosome, is
> resolved to DEL/DUP/INS by `find_edited_reads.py:420` and never reaches the BND queue at all — so
> the two "deletion" rows were unreachable-by-construction wrong. Confirmed two ways: all eight
> events carry both junctions of the reciprocal pair (a deletion makes one), and there are zero
> collinear junctions at any of the eight cut pairs even unfiltered. The shipped labels are now
> `multi-cut inversion` / `inversion at cut site`. Kept here as written for the record.


### Not done

- Not wired into the Nextflow pipeline. It runs standalone over `*.offtarget_analysis.tsv`.
- The cohort numbers above come from **pre-fix** inputs. A full re-run with the corrected caller
  would drop the 4 promiscuity drops (removed upstream instead) and leave the 25-junction queue
  unchanged. Worth confirming when the cohort is next re-run for another reason; not worth a
  dedicated 32-sample run.

---

## 2026-08-10 — Item 1: consolidate rules 1–3 into the caller. **Scoped down, deliberately.**

The plan asked to move review rules 1–3 into `find_edited_reads.py`. Measuring first changed what
was worth doing. Two real defects surfaced; neither is fixed by moving a rule, and the fixes that
*would* move rules turn out to make things worse.

### Finding A — rule 2 is duplicated *inside* the caller, and the two copies disagree

There are three distance tests, not one:

| where | formula | used for |
|---|---|---|
| `:2035`, `:2057`, `:2085`, `:2107` | `distance > max AND distance2 > max` | per-read rejection |
| `:2193` → filtered at `:2229` | `min(\|pos−PAM\|, \|pos+len(ref)−1−PAM\|)` — capital `Distance` | **the site filter** |
| `:593` | `min(\|pos−PAM\|)`, anchor base only — lowercase `distance` | **what is reported** |

`Distance` is a minimum over a superset, so `Distance <= distance` always. The caller therefore
**filters on one quantity and reports a larger one**. Consequence, measured on the cohort audit
trail: **32 rows report `cut_dist_min > 25` under a nominal `-d 25` cap.** Downstream cannot
reproduce the caller's own rule-2 decision from the reported field.

### Finding B — `-x` discards evidence without recording it

`-x/--max-in-control` (default 0) removes control-supported events before they are written, so a
site reports `indel_count = 5` whether it suppressed one event or fifty. Measured on the 92-site
KLF12 subset:

| | `-x 0` (default) | `-x` disabled |
|---|---|---|
| `indel_count` | 90 | **170** |
| `indel_reads` | 412 | **2,966** |
| `bnd_count` | 40 | 41 |

**81 events and 2,554 indel reads suppressed across 42 of 92 sites, with no record.** That is the
caller performing most of the germline removal *before* `review_filter`'s rule 1 ever runs.

Confirmed at the same time: `control_indel_reads` and `control_bnd_reads` are **identical with `-x`
on or off**, so the summation-before-filter fix works and the plan's second acceptance criterion —
"control columns always populated regardless of `-x`" — was already satisfied.

### What was implemented: two audit columns, no behavioural change

| column | meaning |
|---|---|
| `n_control_filtered` | events at this site removed by `-x` |
| `min_cut_distance` | the site minimum of `Distance` — the quantity actually filtered on |

Header is now 19 columns; both writers verified at 19. Purely additive.

**Verification, 92-site KLF12 A/B on identical CRAMs:**

```
cols 17 -> 19 | added ['n_control_filtered', 'min_cut_distance']
pre-existing columns differing : NONE - bit-identical
CROSS-CHECK  n_control_filtered 81  vs  suppressed 81  ->  AGREE
sites where -x suppressed something : 42 / 92
min_cut_distance where events exist : n=68, min 0, max 25
```

`min_cut_distance` maxes at exactly 25, matching `-d`, which is the point: the filtered quantity
never exceeds the cap, while the *reported* per-event `distance` does on 32 cohort rows. `-1` is the
sentinel for a site with no events. Test suite `30 passed`.

### What was deliberately NOT done, and why

**Rule 2 was not consolidated into the caller.** Tightening `-d` from 25 to 10 would delete the
audit trail's single largest drop category: of 1,498 gated cohort rows, **753 sit at 11–25 bp** and
680 are dropped by `review_filter` as "far from PAM". Filtering them in the caller means they are
never written and the reviewer can never see what was removed. The coarse-then-fine arrangement is
correct; the defect was that the two stages measured different things, which `min_cut_distance` now
exposes without changing behaviour.

**Rule 3 was not moved.** It requires ≥3 distinct indel lengths, which a site with 2 indel reads can
never satisfy — it is a read-depth threshold in disguise (see the threshold-change entry: 837 rows,
zero survivors). `review_filter` computes it from `indel_info` at cohort level where the depth
context is visible. Moving it into the caller adds depth-sensitivity and buys nothing.

**Rules 4 and 5 cannot move**, as the plan already noted: the panel of normals needs every unedited
sample in one pass, and `find_edited_reads.py` runs per-sample.

**`-x` was left at its default.** Flipping it to report-only is the change that would truly
"let downstream decide", but it multiplies `indel_reads` by 7 and would require re-validating the
whole review filter and the panel of normals. It belongs with the off-by-one indel question as a
separate, measured piece of work — not bundled here.

### A bug the cross-check caught, worth recording

The first build of these two columns shipped **empty** — `n_control_filtered` 0 at every site,
`min_cut_distance` -1 at every site — while every pre-existing column stayed bit-identical. The A/B
would have passed on the "nothing regressed" criterion alone.

Cause: the values are computed just before the `-x` filter, but the "Process indel results" block
below it re-initialises its locals, and the initialiser added there ran *after* the computation and
blanked both. `control_bnd_reads` survived the same pattern only because it is reassigned further
down.

It was caught by an independent cross-check, not by the diff:

```
n_control_filtered  =  0     (what the new column claimed)
(-x disabled) minus (-x 0)  =  81     (what was actually suppressed)
```

**Lesson for additive columns: "pre-existing columns unchanged" is necessary but not sufficient.**
A new column also has to be checked against something computed a different way, or an empty column
looks exactly like a clean result — the same failure mode as the PoN coverage guard in
`review_filter.py`, and the same one that made `control_bnd_reads` read zero for a whole afternoon.

---

## 2026-08-10 — Item 4: DRAGEN systematic-noise panels

### Part A — SV panels for breakends

Tested all three panels against the 36 gated junctions from item 3 (25 real queue, 11 artifacts).
The metric is **discrimination**, not coverage: a panel is only useful if it flags artifacts at a
higher rate than real junctions.

| panel | records | real queue flagged | artifacts flagged |
|---|---|---|---|
| `WGS_hg38_v3.1.0` | 311,395 | **0 / 25** | **11 / 11** |
| `IDPF_WGS_v3.0.0` | 2,626,364 | 25 / 25 | 11 / 11 |
| `WGS_FF_Heme_v3.1.0` | 2,195,842 | 24 / 25 | 11 / 11 |

**IDPF and FF_Heme are not merely useless — they are harmful.** Across the full 2,403-junction set
they flag the real queue at a *higher* rate than the noise (IDPF: 100% of real at every slop and
both matching modes). Used as a blacklist they would preferentially delete findings. The plan
already warned about IDPF at 200 bp; the measurement is worse than that — it fails at 0 bp too.

**`WGS_hg38_v3.1.0` is perfect on this cohort**, and independently: it flags **every one of the 11
artifacts and none of the 25 real junctions**, using no PAM distance, no cohort context and no
matched control. That makes it the only rule here that works on a **single sample on day one**,
which is exactly the gap the promiscuity rule cannot cover.

Slop is load-bearing:

| slop | real flagged | artifacts flagged |
|---|---|---|
| 0, either end | 0 | 11 |
| 50, either end | 0 | 11 |
| **200, either end** | **3** | 11 |

Default set to **50**; the help text says not to raise it to 200.

**Wired into `review_filter_bnd.py`** as optional `--sv-noise` / `--sv-noise-slop`, on-target
exempt. Verified: with and without the panel the queue is **identical (25 = 25)**. It reports
`dropped, DRAGEN systematic noise: 0` because `np.select` is first-match-wins and all 11 are already
attributed to rules 2 and 3 — the panel is **corroborating, not additive**, on this data. Its value
is (a) independent confirmation that every drop is a known artifact locus, and (b) a rule that still
works when there is one sample and no cohort.

```bash
review_filter_bnd.py <tsvs> --sv-noise \
  /storage2/fs1/dspencer/Active/shared/refdata/hg38/dragenfiles/WGS_hg38_v3.1.0_systematic_noise.sv.bedpe.gz \
  -o bnd_queue.tsv
```

### Part B — the SNV noise BED

`IDPF_WGS_hg38_v.2.0.0_systematic_noise.snv.bed.gz`, 1 GB, BGZF (tabix-able but shipped without an
index). Format is 7 columns: `chrom start end mean max alleles n_samples`, and the allele column
carries `D` codes alongside `A/C/G/T`, so it is **not purely an SNV panel** — worth testing after
all, contrary to the plan's prior. Provenance matters: `##PON SAMPLES` lists ~46 leukemia normals,
so it is a heme panel applied to CAR-T (T-cell) data.

Streamed the whole file once (97,851,753 records) against all 1,498 gated indel rows, ±2 bp.

**Allele vocabulary settles the first question:** `G` 30.3M, `C` 29.8M, `T` 15.1M, `A` 14.8M,
**`D` 1,299,670**, **`I` 892,925**, plus multi-allelic combinations. The `D`/`I` codes mean this is
**not purely an SNV panel** — it carries deletion and insertion noise, so the plan's prior ("expect
a negative result, it is an SNV map not an indel map") was wrong and the test was worth running.

| class | rows | flagged |
|---|---|---|
| dropped (artifacts) | 1,410 | **467 (33.1%)** |
| off-target keeps | 7 | 1 (14.3%) |
| on-target keeps | 81 | **1 (1.2%)** |

**27x enrichment for artifacts over real on-target edits.** That is genuine discrimination — far
better than the IDPF *SV* panel managed, and comparable in shape to the repeat masks already used by
rule 5.

**It fails the strict acceptance bar, but only just:** it flags 1 of 81 on-target sites,
`IKZF2-KO-DNA chr2:213,147,790` — which is the IKZF2 cut site itself, and also one end of that
sample's multi-cut deletion. A real cut site can sit in a locus a generic panel calls noisy, for the
same mappability reasons that made it a predicted off-target in the first place. So this panel
**cannot be used as an unconditional blacklist**; it needs the same `is_off` exemption that rules 4
and 5 already apply, and with that exemption it touches no on-target site by construction.

### The two flagged keeps are not equivalent — and one is a real challenge to a finding

| site | class | noise mean | noise max | alleles | PoN samples |
|---|---|---|---|---|---|
| `CART_NS0027-B2M_2` chr1:28,580,333 | **off-target keep** | 0.0136 | **0.0784** | **`C,D`** | **17** |
| `IKZF2-KO-DNA` chr2:213,147,790 | on-target keep | 0.0006 | 0.0269 | `G` | 1 |

The on-target hit is a coincidence and dismissible on its own fields: **one** PoN sample, SNV-only
(`G`), max 2.7% — against an edit called at **74% VAF**. Nothing to answer.

The off-target hit is substantive. `chr1:28,580,332` is one of the **7 off-target candidates in the
deck**, one of the two low-VAF "hard cases", called at **9.4% VAF**. The panel says that locus is
recurrently noisy in **17 of ~46 leukemia normals**, that the noise includes a **deletion (`D`)**,
and that it reaches **7.8%** there. Our candidate sits barely above that. This is independent
evidence — from a panel built on different donors, different chemistry, different lab — that the
B2M chr1 candidate may be systematic noise rather than an edit. **It does not settle it, but it is
the strongest external challenge to any of the 7.**

Note this is the same site the earlier `pon_min_reads` work flagged as marginal: control VAF 0.0085,
1 control read, and it would have been blacklisted at `pon_min_reads = 1`.

**Rule design follows from those two rows.** A bare interval hit is the wrong test — it fires on the
IKZF2 coincidence too. Sweeping both candidate conditions against the queue shows which one is
actually load-bearing:

| `min_donors` (with `D`/`I` required) | dropped rows hit | off-target keeps hit | **on-target keeps hit** |
|---|---|---|---|
| 1 | 53 | 1 | **0** |
| 3 | 12 | 1 | **0** |
| 5 | 10 | 1 | **0** |
| 20 | 1 | 0 | **0** |

**The `D`/`I` allele requirement alone excludes every on-target site, at every threshold including
1.** IKZF2's alleles are `G` — SNV-only noise, which says nothing about an indel call at the same
coordinate. Only 25 of the 359 panel records touching the query set carry `D`/`I` at all, so the
allele test is doing the protective work; the donor floor is a confidence knob controlling how many
already-dropped rows are additionally corroborated (53 → 12 → 10 going 1 → 3 → 5).

Default set to **`--snv-noise-min-donors 3`**.

**Recommendation.** Usable as a *supplement* to the run's own PoN for indels, never a replacement:
the self-built panel is guide-matched and donor-matched, this one is 46 leukemia normals. Worth
adding to `review_filter.py` as an optional `--snv-noise` rule behind the on-target exemption, and
worth measuring against the queue before enabling by default — it flags 1 of the 7 off-target
candidates, so it would change findings, not just drop noise. **Implemented** as optional rule 6 in `review_filter.py`
(`--snv-noise`, `--snv-noise-min-donors`), off by default. It is placed **last** in the `np.select`
order so first-match-wins guarantees it cannot re-attribute any existing drop — it can only claim
rows the five existing rules kept. Verified: without the flag the run is unchanged (88 queue, same
per-rule counts, `systematic noise: 0`).

### Rule 6 in production — measured effect

```
review_filter.py <tsvs> --pon <pon> --repeats <beds> --min-reads 2 --min-vaf 0.005 \
                 --snv-noise IDPF_WGS_hg38_v.2.0.0_systematic_noise.snv.bed.gz
```

| | without rule 6 | with rule 6 |
|---|---|---|
| queue | 88 | **87** |
| on-target | 81 | **81** (unchanged) |
| off-target | 7 | **6** |
| rows added | — | none |

The single row it removes is **`CART_NS0027-B2M_2 chr1:28,580,332`**, dropped as
`systematic noise (external panel)`. Every on-target site is retained, and nothing is added or
re-attributed.

**This changes a finding, so it is off by default.** Enabling it takes the off-target list from 7 to
6, and the site it removes is one of the two low-VAF hard cases. The evidence for removing it is
strong — 17 of ~46 unrelated donors show deletion noise at that exact locus, reaching 7.8% against
our 9.4% call — but it is an external heme panel judging a CAR-T sample, so the call belongs to
whoever is signing off on the cohort, not to a default.

Cost: streams the 1 GB panel once per invocation, ~3 minutes. The file is BGZF and the directory is
writable, so `tabix -p bed` would make this instant; not done, to avoid writing into shared refdata
without asking.

---

## Wiring into the pipeline (2026-08-11)

Both rule 6 and the breakend filter were implemented as Python but were **not reachable from the
pipeline**: `modules/local/review_filter.nf` never passed `--snv-noise`, and `review_filter_bnd.py`
had no module, no subworkflow call and no params. A cohort rerun would have silently reproduced
`results_cart_ponfix`. That gap is now closed.

### What changed

| file | change |
|---|---|
| `nextflow.config` | `review_snv_noise`, `review_snv_noise_min_donors`; the `review_filter_bnd` / `review_bnd_*` / `review_sv_noise*` block |
| `nextflow_schema.json` | all of the above registered under `review_filter_options`; stale defaults corrected (`review_min_reads` 10→2, `review_min_vaf` 0.05→0.005, `pon_min_reads` 3→2) |
| `modules/local/review_filter.nf` | fourth input `path snv_noise`, threaded into **both** script invocations |
| `modules/local/review_filter_bnd.nf` | new process `REVIEW_FILTER_BND` |
| `conf/modules.config` | `REVIEW_FILTER_BND` publishes to `${params.outdir}/review` |
| `subworkflows/local/scge_analysis.nf` | `ch_snv_noise`, `ch_sv_noise`, and `ch_analysis_tsvs` shared by both processes |

`ch_analysis_tsvs` is a `.collect()`, which yields a *value* channel, so the same staged list feeds
both review processes without being consumed by the first.

The schema defaults were genuinely wrong before this, not merely untidy: nf-core lint compares
`nextflow_schema.json` defaults against `nextflow.config`, and all three had drifted when the gate
was lowered.

### Verified before committing

Both wired command lines were run against the 32 `results_cart_ponfix` analysis tables using the
exact arguments the modules emit.

* `REVIEW_FILTER` with rule 6 → **87** rows (baseline 88), `is_target` 81 on-target / 6 off-target.
  The one removed row is `chr1:28,580,333`; nothing added. Note this is the 1-based `end`; earlier
  notes call the same site `28,580,332`, which is its 0-based start.
* `REVIEW_FILTER_BND` → **25** rows: 2,404 junctions → 36 gated → 25 queued (23 multi-cut deletions
  in 7 samples spanning 413–125,747 bp, 2 deletions at a cut site). `dropped, DRAGEN systematic
  noise: 0` — the WGS v3.1.0 panel flags none of the real junctions, which is the expected and
  desired behaviour for the one safe panel.
* `nextflow run . -preview` with both panels wired builds the DAG with `REVIEW_FILTER` and
  `REVIEW_FILTER_BND` present and completes successfully.
* `nextflow run . --help` lists every new param — the real test that schema registration is correct.
* `pytest tests/` → 30 passed.

`-profile stub -stub-run` cannot be used as a wiring check here: the stub profile uses the local
executor without containers, so `PARSE_INPUT_SAMPLESHEET` (which has no stub block) dies on a
missing pandas long before the review processes. This is pre-existing and unrelated. `-preview` is
the check that works.

### Defaults chosen, and why

`review_filter_bnd` defaults to **true** — it is cheap and additive, producing a second queue file
without touching the indel result. `review_snv_noise` and `review_sv_noise` both default to **null**,
so the rules stay off unless a panel is passed. Rule 6 changes a finding (7 off-target → 6), and the
SV panel choice is dangerous enough to require a deliberate act: `IDPF_WGS_v3.0.0` flags 25/25 real
junctions and `FF_Heme_v3.1.0` flags 24/25, so wiring either would erase the entire breakend result.
Only `WGS_hg38_v3.1.0` is safe. If the breakend queue ever collapses toward zero, suspect the BEDPE.

### Tabix: deliberately not done

`load_snv_noise` streams the file line by line, so an index buys nothing without also adding a
pysam/tabix code path. Against a multi-hour cohort run the ~3 minutes (twice, since the module runs
the script for both the queue and the audit) is not worth a new code path plus writing into shared
refdata. Recorded as a known optimization, not a defect.
