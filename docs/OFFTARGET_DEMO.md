# Off-target CRISPR edit finder — a 5-minute demo

A hands-on tour you can run on a login node in about a minute. Nothing here needs cohort
data, a cluster allocation, or the DRAGEN reference: every command below runs against a
**5.5 KB synthetic genome** committed to the repo.

For the full scientific writeup — what the workflow is and is not, how the shape model was
trained, what the AAVS1 cohort actually showed — see
[`OFFTARGET_WORKFLOW.md`](OFFTARGET_WORKFLOW.md). This page is the "show me it working" one.

---

## What the finder does, in three sentences

CRISPR cuts where you aimed it, and sometimes where you didn't. This workflow takes deep
error-corrected sequencing (ECS) at a panel of candidate sites and finds the reads that
carry an edit, then asks whether ordinary 30x whole-genome sequencing of the *same* sample
could have found those same edits on its own.

Two arms, one join: **ECS is the truth arm**, **WGS is the arm under test**.

---

## Setup (once)

```bash
# nextflow + apptainer on the RIS cluster
S=/storage2/fs1/dspencer/Active/spencerlab/apps/modules/spack/software/linux-rhel9-x86_64_v4/gcc-13.3.0
export PATH="$S/nextflow-25.10.4-fz3mazhspz5uruz5ymyxhcwcauohjgfm/bin:$S/apptainer-1.4.5-ifu7outvqr62l343tqfxuxcguzftiuqe/bin:$PATH"

# nf-test (once, anywhere on PATH)
curl -fsSL https://code.askimed.com/install/nf-test | bash
```

---

## Demo 1 — the whole test suite, real execution, ~1 minute

```bash
nf-test test --profile stub,apptainer
```

```
OFFTARGET_METRICS
  Test 'computes the real metrics and emits every declared channel'            PASSED
  Test '--betas reaches the script: a custom beta set changes the columns'     PASSED
  Test '--negatives all_label0 is honoured as a distinct negative set'         PASSED
ECS_INDELS
  Test 'calls the edit and reports the right VAF — tagged BAM off by default'  PASSED
  Test 'offtarget_tagged_bam=true emits a per-read tagged BAM and its index'   PASSED

SUCCESS: Executed 5 tests
```

These are not stubs. Each one starts the real container and runs the real Python.

> **Why that distinction matters.** The repo also has `-stub-run` tests, and they are
> useful — but a stub process body only `touch`es its output files. It never executes a
> line of `bin/*.py`. A module that passes a flag the script does not accept stays green
> through every stub test in the repo. That exact bug class once survived CI and then died
> three hours into a cluster run. nf-test is the tier that catches it.

---

## The toy genome

Everything below runs on `tests/fixtures/ecs/` — 5.5 KB, regenerable with
`python tests/fixtures/make_ecs_fixture.py`:

```
chr1, 3000 bp
  target A @ 1001   the EDITED site  — 12 of 20 evaluable read pairs carry a 5 bp deletion
  target B @ 1201   a QUIET site     — no edit

edited.cram    35 pairs / 70 records
control.cram   20 clean pairs (the matched normal)
```

The 35 pairs are awkward on purpose:

| pairs | what they are | what should happen |
|---|---|---|
| 12 | carry a 5 bp deletion | counted as edited |
| 8 | clean, span the target | counted as reference |
| 2 | duplicate-flagged | **excluded** from the denominator |
| 2 | MAPQ 3 | **excluded** (floor is 20) |
| 2 | NM 6 | **excluded** (ceiling is 4) |
| 9 | sit where the two ±150 bp windows **overlap** | visited by both targets, must be counted once |

That last row is the interesting one: targets A and B are 200 bp apart, so their fetch
windows overlap. Any per-read output has to survive being visited twice.

---

## Demo 2 — find the edit

```bash
IMG=/storage2/fs1/dspencer/Active/spencerlab/abonney/apptainer_cache/ghcr.io-dhslab-docker-scge-offtarget-260710.img
F=tests/fixtures/ecs

apptainer exec -B /storage2 "$IMG" python bin/find_edited_reads.py \
    --fasta $F/ref.fa \
    --edited-bam $F/edited.cram --control-bam $F/control.cram \
    --target-file $F/targets.vcf \
    -o demo.tsv
```

```
chrom  start  end   total_reads  indel_reads  indel_fraction  control_reads  control_indel_reads
chr1   1000   1001  20           12           0.6             20             0
chr1   1200   1201  23           0            0.0             20             0
```

Read that as: **at the edited site, 12 of 20 evaluable reads carry the deletion (VAF 0.60),
and the matched control has none of them.** That last column is what makes the call somatic
rather than germline. The quiet site stays quiet.

Note `total_reads` is 20, not 26 — the duplicate, low-MAPQ and high-mismatch pairs were
dropped before the denominator was formed. If a filter ever regresses, that number climbs
and the VAF silently falls. The test asserts on it for exactly that reason.

---

## Demo 3 — per-read tags for IGV *(new)*

Add one flag and every read gets an `XC` tag naming how the caller classified it:

```bash
apptainer exec -B /storage2 "$IMG" python bin/find_edited_reads.py \
    --fasta $F/ref.fa \
    --edited-bam $F/edited.cram --control-bam $F/control.cram \
    --target-file $F/targets.vcf \
    --tagged-bam-out demo.tagged.bam -o demo.tsv
```

```
  31  Unedited_WT
  15  Skipped_NoSpan
  12  Edited_Deletion_5bp
   4  Skipped_Duplicate
   4  Skipped_LowMapQ
   4  Skipped_Mismatches
──────────────────────────
  70  records, 35 unique read names
```

Load `demo.tagged.bam` in IGV and use **Color alignments by → tag → XC**. A reviewer can now
see *why* a read was or wasn't counted, instead of taking the caller's word for it.

Two things worth pointing at on a slide:

- **70 records, 35 unique names.** Every read is written exactly once, even the 9 pairs in
  the window overlap that two different targets both visited.
- **The `Skipped_` classes are visible.** Reads that were excluded are in the file and
  labelled, not silently missing.

In the pipeline this is `--offtarget_tagged_bam`, and it is **off by default**: on a real
AAVS1 sample (1149 targets, ~11000x) it produces an 800 MB BAM and pushes the caller's peak
memory from 1.4 GB to 6.3 GB. Turn it on for review, not for a cohort sweep.

---

## Demo 4 — clinical metrics *(new)*

```bash
apptainer exec -B /storage2 "$IMG" python bin/offtarget_metrics.py \
    --training tests/fixtures/training_mini.tsv \
    --betas 2,5 --out-json m.json --out-txt m.txt && cat m.txt
```

```
-- RANKING (score as a continuous ranker; ECS-label denominator) --
  n_pos / n_neg      : 6 / 4   prevalence 0.6000
  PR-AUC             : 0.8552   (1.43x prevalence)
  excluded, no score : 1 (1 pos / 0 neg) — INSUFFICIENT COVERAGE, never filled with 0

-- OPERATING POINT (verdict contains 'LIKELY EDIT'; ECS-label denominator) --
  TP 5   FP 2   FN 1   TN 2
  precision : 0.7143      recall : 0.8333      recall incl. unevaluable : 0.7143
  F1 0.7692    F2 0.8065    F5 0.8280
```

Three things this output is built to stop you getting wrong:

**1. F-beta is recall-weighted, and you can see it.** Here recall (0.833) beats precision
(0.714), so F rises with beta: F1 < F2 < F5. For a screening assay a missed off-target edit
costs far more than a followed-up false one, which is the whole argument for reporting F2/F5
rather than F1. F-beta is monotone in beta, so F1 is always an *endpoint* — if you ever see
F1 in the middle, the weights got applied backwards.

**2. The uncovered site is excluded, not counted as zero.** One credible positive had no WGS
coverage at all. It is dropped from the ranking and reported separately, and `recall incl.
unevaluable` (0.7143) is printed next to `recall` (0.8333) so the depth floor can never hide
inside a good-looking number.

**3. Precision is against the ECS label, never against manual review.** The report says so
out loud, in the file. Manual review has 55 confirmed positives and 80,440 NaNs, where NaN
means *unreviewed*, not *rejected* — it contains no confirmed negatives. Any precision
computed from it would score every discovery the reviewers never reached as a false
positive. Recall against manual review is reported separately, and as recall only, by
`bin/validate_recall.py`.

---

## Where the tests live and what each tier is for

```bash
bash run_offtarget_tests.sh        # all five tiers, ~3 min
```

| tier | what it runs | what it proves | catches |
|---|---|---|---|
| 0 | `ast.parse`, `nextflow -preview` | it is syntactically valid | typos |
| 1 | `pytest tests/*.py` | **the science is right** | wrong math, wrong joins |
| 2 | unpickle the shape model | the sklearn pin holds | version drift |
| 3 | `nextflow -stub-run` | **the DAG wires up** | broken channels |
| 4 | `nf-test` | **the two halves connect** | module ↔ script contract |
| 5 | real AAVS1 cohort (SLURM) | it works on real data | everything else |

Tiers 1 and 3 are the two halves that tier 4 joins: pytest runs the Python but never through
Nextflow; the stub run goes through Nextflow but never runs the Python.

Run one feature's tests only:

```bash
nf-test test --tag metrics     --profile stub,apptainer
nf-test test --tag tagged_bam  --profile stub,apptainer
```

---

## Honest limits

- The fixtures are **synthetic**. They prove the plumbing and the arithmetic, not that the
  caller is right about real CRISPR biology — that is what the AAVS1 cohort run is for.
- The metrics fixture is hand-built to produce a specific confusion matrix. Its PR-AUC is
  not a performance claim about the model; it is a fixed number chosen so the test can
  assert on it.
- Two blind spots remain unmeasured on real data: the **sub-5% VAF floor**, and a **de-novo
  off-target positive control**. Neither is addressed by anything on this page.
