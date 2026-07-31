# Off-target edit finder — what it does, how it works, what it does not yet prove

Written for the questions that actually get asked. Every number here is reproducible from the
commands at the bottom; nothing is quoted from memory.

---

## 1. What features decide "CRISPR edit" vs "sequencing error"?

Seven, computed per candidate site in `bin/features.py:184`:

| feature | what it measures | why it separates edits from noise |
|---|---|---|
| `indel_frac` | fraction of spanning reads carrying an indel | an edit is present in a real fraction of the cells |
| `conc_ratio` | fraction of indel reads within ±2 bp of the modal position | a cut happens at **one** place |
| `pos_conc` | same, as a fraction of *all* spanning reads | clonality relative to total depth |
| `pos_mad` | median absolute deviation of indel positions | noise scatters; cuts do not |
| `modal_len` | the dominant indel length | |
| `modal_mapq` | median mapping quality of the indel reads | repeats mismap and look like edits |
| `softclip_frac` | soft-clipping at the site | the signature of a structural break |

Three more were added after auditing the model against the manual-review rules
(`features.py:MODEL_FEATURES` is now 10 features):

| feature | what it measures | the review rule it encodes |
|---|---|---|
| `cut_dist` | bp from the observed indel to the **predicted cut site** | "is the edit near the PAM?" |
| `homopolymer_len` | longest single-base run touching the site | "is it in a repeat?" |
| `repeat_frac` | fraction of a ±25 bp window inside a 1–4mer tandem repeat | "is it in a repeat?" |

`cut_dist` is NaN, not 0, when no indel was observed — 0 would falsely claim the indel sits
exactly on the cut. `HistGradientBoostingClassifier` handles NaN natively.

**The one-sentence version:** a real edit is *many reads agreeing on one position, at the cut,
outside a repeat*; noise is *a few reads disagreeing everywhere, or agreeing inside a repeat*.

### What the original 7 features were missing

Being straight about this is better than being caught by it:

1. **Distance to the PAM / cut site — was missing entirely, now a model feature.** Cas9 cuts
   ~3 bp from the PAM, so a genuine edit sits within a few bp of the predicted site. Validated
   **held-out** on 146 human-confirmed CART edit rows (the feature was never in training): median
   distance 1 bp, 98.6% within 20 bp, while sites called at curated-negative AAVS1 loci run to a
   median of 24 bp. See §3.

   Note it is now a *feature*, not the hard gate first shipped. A gate destroyed 2 genuine CART
   edits (`PLCB2 chr15:40302475` at 108 bp with `indel_frac` 0.613, and `CTLA4 chr2:203870831` at
   58 bp) because a threshold cannot trade distance off against strong evidence. A feature can.
2. **Indel length heterogeneity — now measurable, but not yet a validated discriminator.** Your
   review rule is that *if every read shows the exact same indel length, it is probably a
   sequencing/alignment artifact*: real CRISPR pileups carry a mix of lengths at slightly
   different positions, all overlapping one cut site. The model has `modal_len` (the mode) and
   nothing about the *spread*; the length array is computed and discarded in `features.py`.

   The subtlety is the interesting biology: this is *not* `conc_ratio`.
   **Position should be concordant; length should not.**

   It is now measurable directly from the tagged BAMs (`bin/indel_length_diversity.py`), and the
   AAVS1 on-target sites look exactly as you predicted — 19 and 16 *different* indel lengths at
   one cut position (1 bp, 3 bp, 10 bp, 12 bp, 13 bp, 14 bp, 17 bp, 33 bp, 38 bp ...). That is a
   good thing to show live in IGV.

   **But it does not yet separate, and the raw number must not go on a slide unqualified:**

   | | TP (n=2) | TN (n=26) |
   |---|---|---|
   | raw distinct lengths | 19, 16 | median 2, **max 14** |
   | depth-controlled at 30 edit reads | 8, 9 | **9**, 5, 2 |

   The raw count is confounded with depth — the real edits carry 119–149 edit reads, most
   negatives carry 1–17. Rarefy every site to the same 30 reads and the best negative
   (`chr6:36797655`) scores 9, matching the two real edits. Adding it as a model feature would
   need a retrain **and** depth conditioning; the right form is the raw count fed alongside
   `n_edit_reads` so the tree model learns the interaction — *not* normalised Shannon entropy,
   which divides out the signal (TP normalised entropy is 0.55–0.59 while several negatives reach
   1.000, inverting the ranking).

3. **Repeat context — was only an indirect proxy, now direct.** `modal_mapq` catches repeats
   that *mismap*, but short homopolymers and STRs still map uniquely at high MAPQ, so they were
   invisible. `homopolymer_len` and `repeat_frac` read the reference directly and need no external
   annotation track. This encodes the rule you check first in IGV.

---

## 2. What kind of model is it, and how does it work?

`assets/models/wgs_shape_model.pkl`:

- **`HistGradientBoostingClassifier`** (scikit-learn), `max_depth=3` — a small gradient-boosted
  tree ensemble, not a neural network.
- **`role: stage2_shape_ranker`.** Stage 1 is a cheap rule filter; stage 2 is this model ranking
  what survives.
- **Trained on 51 positive and 383 negative loci** from the CART ECS cohort.
- **Positive-unlabeled**: unreviewed sites were treated as artifacts, so some "negatives" are
  really just unreviewed.
- **Depth-augmented**: training reads were subsampled to simulate 0/20/30/50x, so the model
  behaves at 30x WGS rather than only at ECS depth.

Two weaknesses worth volunteering before someone finds them:

- **51 positives is a small training set.** Treat it as a prioritiser, not an oracle.
- **It is a ranker, not a detector.** It demotes some unambiguous edits. That is exactly why the
  **high-evidence rescue** exists (`bin/score.py:145-147`): any site with `indel_frac >= 0.15`,
  `conc_ratio >= 0.5` and `spanning >= 20` is called regardless of model score. The rescue runs
  *after* the germline and low-MAPQ vetoes, so it can never resurrect a germline variant.

---

## 3. Does the ML actually beat just running the SCGE indel finder?

This is the honest core of it. The SCGE pipeline's `get_indels.nf` already finds indels — it runs
`bin/find_edited_reads.py` with the ML commented out. Its cost is the **size of the review queue**
a human then has to open in IGV.

Scored against the curated AAVS1 truth (`crispr_ml/AAVS1_training_{tp,tn}.tsv`: 2 confirmed edit
sites, 6,035 confirmed non-edit sites), collapsed to unique sites:

| arm | queue | TP | FP | recall | precision |
|---|---|---|---|---|---|
| ECS rules, any indel read | 4,636 | 2 | 4,602 | 1.000 | 0.000 |
| ECS rules + somatic gate | 398 | 2 | 381 | 1.000 | 0.005 |
| **WGS rules only, no model** | **26** | 2 | 19 | 1.000 | **0.095** |
| WGS + ML shape ranker | 35 | 2 | 28 | 1.000 | 0.067 |
| **WGS + ML + `cut_dist<=10`** | **10** | 2 | 5 | 1.000 | **0.286** |

**What this shows: a 40x smaller review queue than the current workflow, at identical recall.**
398 sites to open in IGV becomes 10.

**What it does not show — and you should say so first:** on *this* cohort the ML does **not** beat
simple rules. Rules-only gives 26 sites at precision 0.095; adding the model gives 35 at 0.067.
The large win comes from the cut-site gate, which is a **rule**, not the model.

### Why this cohort cannot settle the question

The model exists to find **low-VAF** edits that rules threshold away. This truth set cannot test
that:

- Both curated positives are **on-target**, at WGS `indel_frac` 0.85 and 0.89 — roughly 17x above
  the rule threshold. Any method finds them.
- **Zero** curated positives sit below `indel_frac` 0.05, which is the entire regime the model was
  built for.

The model *is* reaching into that regime: 24 calls that rules would reject, 12 with independent ECS
support at 0.6–3.8% VAF. But when those 12 are checked against the cut-site heuristic, **only 1 of
12 is within 10 bp of the predicted cut and 9 of 12 are beyond 20 bp** — so most are probably
artifacts, and their ECS support is itself in the ECS assay's noise band. Resolving them needs
manual review, which is what the IGV session is for.

**Bottom line to state plainly:** recall and review-burden are demonstrated; the model's low-VAF
advantage is *not yet* demonstrated, because no curated low-VAF positives exist. That is a gap in
the truth set, not a result against the model.

---

## 4. Why is `offtarget_metrics.txt` reporting recall 0.033?

Because that file scores against the **ECS label**, not against curation, and that denominator is
contaminated:

- 3,739 of ~4,600 `ecs_is_edit==1` rows are below 0.5% VAF — ECS assay noise.
- Of the 116 above 5% VAF, **94 carry the same indel in the matched unedited normal** (`ctrl_if` up
  to 1.0, many at `ecs_if` = 1.0000, i.e. homozygous). Those are **germline variants**.

The WGS arm rejects them as `GERMLINE/ARTIFACT (in normal)` — correctly. Counting them as misses
would mean demanding the pipeline report germline variants as CRISPR edits; "fixing" that number
would make the pipeline worse. `bin/validate_recall.py`'s own docstring makes the same point about
the CART cohort.

**Recall against curation is 2/2 = 1.000. Recall against the raw ECS label is not a meaningful
number.** Use `bin/validate_recall_aavs1.py`.

---

## 5. How would a clinician actually use this?

They have a guide RNA and want assurance there are no off-target edits.

```bash
# 1. guide RNA -> predicted off-target sites (Cas-OFFinder), no genome-wide compute
nextflow run . -entry HOTSPOTS --input grna.csv --outdir hotspots

# 2. score the sample at exactly those sites
nextflow run . -entry OFFTARGET --input samplesheet.csv --outdir results
```

Step 1 is the auto-hotspot feature: hand it a spacer sequence and it generates the target list, so
nobody has to build one by hand. Step 2 returns a short ranked list plus, with
`--offtarget_wgs_tagged_bam`, a BAM whose reads are individually labelled so the clinician's
analyst can confirm each call in IGV rather than trusting a score.

For the routine case that is the whole workflow. The genome-wide arm exists for the harder
question — edits at sites with no sequence similarity to the guide (vector integration,
translocations) — which homology prediction cannot enumerate by construction.

---

## 6. Couldn't CRISPResso2 do this?

For a predicted hotspot list, largely yes — and that should be conceded rather than argued.
The differences that matter here:

- It does not use the **25-donor panel of normals**, so germline indels mimicking edits are a
  manual problem. §4 shows that is the dominant contaminant: 94 of 116 high-VAF ECS "edits" are
  germline.
- It does not integrate the orthogonal **ECS 5000x truth labels**, so there is nothing to train or
  validate a ranker against.
- It leans on strict PAM/homology windows, which cannot represent a homology-blind event.

The real answer is the §3 table: same recall, review queue 398 → 10.

---

## 7. Reproducing every number here

```bash
# recall / specificity / precision vs curated truth  (2/2, 0.9940, 0.067)
python bin/validate_recall_aavs1.py \
    --training results_offtarget_aavs1_igv/offtarget/training.tsv \
    --max-cut-dist 10 --require-recall 1.0

# the §3 baseline-vs-ML table
python bin/compare_to_baseline.py \
    --training results_offtarget_aavs1_igv/offtarget/training.tsv \
    --ecs-glob 'results_offtarget_aavs1_igv/offtarget/*.offtarget_analysis.tsv'
```

Both are deterministic and both reproduce identically on the two independent AAVS1 runs
(`results_offtarget_aavs1_rescue` and `results_offtarget_aavs1_igv`).

---

## 8. Known gaps, stated once

- **No curated low-VAF positives**, so the model's core claim is untested (§3).
- **Length heterogeneity is measured but not validated** (§1): the effect is visible and matches
  the biology, but it does not survive a depth control on 2 positives. Needs the CART cohort,
  which has 51.
- **No confirmed off-target edits in AAVS1 at all** — both curated positives are on-target, so
  "off-target recall" is currently unmeasurable in this cohort.
- **Sub-5% VAF floor at 30x WGS is physics**, not tuning: a 1% VAF edit yields ~0.3 supporting
  reads. No model recovers that.
- **De-novo off-target positive control** still missing: the 40 genome-wide `LIKELY EDIT` sites
  that are *not* predicted hotspots are unreviewed. (Of the 49 genome-wide calls, 9 fall on
  predicted hotspots and all 9 are ECS-confirmed — 100% concordance on the overlap.)
- **Base-editor (CBE/ABE) substitutions are not detected.** This pipeline calls indels only.
