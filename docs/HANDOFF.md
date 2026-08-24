# Handoff — off-target arm, `feat/offtarget-wgs`

**Last updated:** 2026-08-21 · **Branch:** `feat/offtarget-wgs` · **Not yet merged to `dev`.**

This is the operational knowledge that was not otherwise written down: how to actually run the
thing, the landmines that cost days, which numbers are real, and what I would do next. The
scientific write-ups live in `docs/OFFTARGET.md`, `docs/NOISE_MODEL*.md` and `docs/PANEL_AS_FILTER.md`
— this file deliberately does not repeat them.

---

## 1. State of play

The pipeline runs end to end on real cohorts. Two arms, both working:

- **Indel arm** — `find_edited_reads.py` → `review_filter.py`. Ran green on all 32 CAR-T WGS
  samples (2026-08-07) and on the 8 AAVS1 samples.
- **Breakend arm** — SA-tag split reads → `review_filter_bnd.py` → `bnd_snapshots.py`. Completed
  2026-08-16, corrected 2026-08-20 (see §5).

Everything is committed and pushed. The working tree is clean.

**The one outstanding action: open the PR.** `feat/offtarget-wgs` is ~50 commits ahead of `dev`
and nothing merges it. <https://github.com/dhslab/nf-core-scge/compare/dev...feat/offtarget-wgs>

---

## 2. How to run it

### 2a. The default workflow (indel + breakend review) — SLURM

This is the non-obvious one. **The review arm lives in `subworkflows/local/scge_analysis.nf`,
which is the DEFAULT workflow (`main.nf:36`), not `-entry OFFTARGET`.** So `run_offtarget.sh` is
the wrong wrapper, and every other `run_*.sh` in the repo is LSF/bsub. The working recipes are
durable at **`logs/run_cart_bnd.sh`** and **`logs/run_cart_nopon.sh.recovered`**:

```
sbatch --partition=condo-dspencer --account=compute2-dspencer ...
  nextflow run main.nf -profile ris2,apptainer \
    -work-dir /scratch2/fs1/dspencer/abonney/cartfull_work -resume \
    --input <samplesheet> --outdir <out> --run_alignment false \
    --offtarget_tagged_bam true --review_repeat_beds ... --review_snv_noise ... --review_sv_noise ...
```

Everything else (`off_target_threshold`, `transgene_*`, `scge_report_qmd`, `run_analysis`) comes
from `nextflow.config` defaults. `params.queue` / `job_group_name` in old params dumps are LSF
leftovers that `ris2` never reads — do not conclude from them that a run was LSF.

**Partition:** `general-cpu` is NOT permitted for this group; use **`condo-dspencer`**. Both the
wrapper *and* the child jobs Nextflow submits need it — pass `--slurm_partition condo-dspencer`
(`nextflow.config:42`, feeds the `ris2` profile's `process.queue` at `:345`). Account `compute2-dspencer`.
Cancelling the parent job does **not** kill the `nf-*` children; `scancel` them explicitly.

### 2b. Environment

- Compute nodes have java 17 at `/usr/bin/java`; **the JupyterLab exec node does not** — that is
  why nextflow tiers SKIP there. There is **no `java17` module**; `module load nextflow/25.10.4`
  brings its own OpenJDK 17.
- Modules: `apptainer/1.4.5`, `nextflow/25.10.4`, `labtools` (this is how you get `Rscript`).
- `export NXF_SINGULARITY_CACHEDIR=/scratch2/fs1/dspencer/apptainer_cache`,
  `export APPTAINER_BINDPATH=/storage2,/scratch2,/home/<user>`.
- **venvs on `/storage2` are not portable across nodes** (`bin/python3 -> /usr/local/bin/python3`
  is per-node). Use the container for anything needing pandas/sklearn.
- `bash script.sh` does not inherit the `module` shell function; re-source `/etc/profile.d/lmod.sh`.
- **`sbatch --output` must point at shared storage, never `/tmp`.** `/tmp` is node-local: the log
  lands on the compute node and is invisible from the submit host — indistinguishable from "the
  job never started". Use `<repo>/logs/<name>.%j.log`. This cost a full misdiagnosis once.

### 2c. Tests

`./run_offtarget_tests.sh` is container-first — pulls `ghcr.io/dhslab/docker-scge-offtarget:260710`
and runs the python tiers inside it (exact prod deps incl. sklearn 1.8.0; pytest injected as a
pure-python overlay since it is not in the prod image). The Nextflow stub tier uses `-profile stub`
(local executor), **not** `-profile ris` (LSF, fails off-scheduler).

`-profile test` **alone fails** (no container engine → `ModuleNotFoundError: pandas`); the
`singularity`/`apptainer` profile is required.

Current state: **nf-test 8/8 green. pytest has 3–6 pre-existing failures** in
`test_offtarget_glue.py` — the missing module is **`pyranges`, not edlib**, plus some `tagged_bam`
cases. They fail with or without any recent change (verified by reverting). Run pytest *inside* the
container to retire them; that has not been done.

---

## 3. Landmines

**Ranked by how much time they cost.**

1. **Never pass `-u` / `--unevaluable-reads` to `find_edited_reads.py`.** It writes a debug log at
   **0.5–1 TB per sample**. It filled `/storage2`, then blew the *shared* `dspencer` group quota on
   `/scratch2` (9.1 TB across 13 concurrent tasks), and killed three CAR-T runs. It is write-only,
   not an emitted output, not consumed downstream. Now gated behind
   `params.offtarget_ecs_unevaluable_log`, **default false**. Leave it false.
   `/scratch2/fs1/dspencer` has a shared group quota (~9–10 TB); `df` shows the whole 1.5 PB
   filesystem and hides it. GPFS accounting lags 1–2 min after you delete.

2. **Editing ANY file in `bin/` invalidates the entire `-resume` cache.** Nextflow stages the whole
   `bin/` directory into every task and folds its hash into the task cache key. A one-line change to
   one script forces a full cohort re-execution — not just the tasks that call it. Verified: five
   changed files in `bin/`, everything else untouched, run reported **cached: 0**. Budget a full
   re-run (~5 h for the 32-sample CAR-T cohort) whenever `bin/` is touched. Do not promise
   "upstream tasks will still cache".

3. **The RIS SLURM controller drops connections, and that kills whole runs.** A failed *submission*
   is fatal to the run, not just to one task — one `Connection reset by peer` lost 160 pending tasks
   after 104 had succeeded. Nextflow's default `executor.retry.reason` matches only
   `Socket timed out on send/recv operation`. The **`ris2` profile widens it** to cover
   "Connection reset by peer" / "Unable to contact slurm controller" / "Connection refused". Copy
   that pattern into any other long grid run here. If `squeue` itself returns
   `Unable to contact slurm controller`, the outage is cluster-side — wait, do not resubmit.

4. **`-stub-run` does not execute `bin/*.py`.** It only touches output files, so it cannot catch a
   Python runtime error. A real run is the only check that does. This is how a `verdict()` arity
   change (2→3 return values) shipped past `py_compile` *and* the stub run and broke the
   genome-wide arm two minutes into a cluster run. Guard added:
   `tests/test_offtarget_glue.py::test_verdict_callers_unpack_three_values` walks the AST of every
   `bin/*.py`.

5. **A Nextflow driver does not exit when the pipeline errors.** It prints "Finishing pending tasks
   before exit" and drains every in-flight task, which can take hours. `squeue` is the truth, not
   the log. Two concurrent runs from the same launch dir share `.nextflow/` and will collide — use
   `sbatch --dependency=afterany:<jobid>` to serialize. Also: running stray `nextflow` commands
   (`config`, `clean`, stub) in a **live launch dir** rotates `.nextflow.log` and can disturb the
   running head.

6. **Verify the artifact, not the exit code.** Five pre-existing bugs on the shared analysis path
   all produced green-looking states with missing or empty outputs (missing `+x` bits → exit 126;
   a hardcoded circos version → exit 127; quarto writing its cache to a read-only `$HOME`;
   `pip install --user` at task time; an Excel writer whose openpyxl fallback was decorative and
   silently skipped every plot). All fixed. Check file counts and embedded-image counts.

7. **`bin/*.py` need the executable bit.** Nextflow puts `bin/` on PATH and invokes scripts bare.
   `git` records 100755 but `core.fileMode` is off here, so `git status` is silent when the bit is
   lost. This has bitten both this repo and `dragenflow`.

---

## 4. Truth sets — which numbers are real

**This is the single most important section.** Most of the "bad results" in this project's history
were truth-set artifacts, not detection failures.

### Never validate against the ECS threshold label

`recall_vs_vaf.py` used to score against `ecs_is_edit`, whose denominator counted **any** nonzero
ECS indel fraction — 47,170 "positives" below 0.5% VAF in the CAR-T run. That produced a reported
recall of ~2% for a pipeline that was actually at 90%. Same trap on AAVS1: the ~3% figure is a
germline-contaminated denominator (of 116 sites above 5% VAF, **94 carry the same indel in the
matched normal**, many homozygous).

**Always validate against manual review.**

### The real truth sources

| cohort | file | notes |
|---|---|---|
| CAR-T ECS | `Manual Indel Review/cart_ecs/cart_ecs_merged.csv.gz` | `manual_review==1`; reaches 2.48% VAF; **positives-only** |
| CAR-T WGS | `Manual Indel Review/cart_wgs/cart_wgs_merged.xlsx`, sheet `gold_wgs` | two-class but ≥5% VAF only |
| AAVS1 | `crispr_ml/AAVS1_training_{tp,tn}.tsv` | **read-level**, not site-level |

**AAVS1 gotchas:** the files are ragged — read with **positional** `usecols=[0,1,2]`; named usecols
silently yields an empty frame. Collapses to **2 TP sites, 6037 TN sites**; both TPs
(`chr19:55,115,732` and `chr19:55,115,752`) are **on-target**. Both TP sites also appear in the TN
set (an edited site contains reference reads too) — **subtract them or you get 2 phantom FPs**.
Join on the table's **`end`** column, not `start` (`start` is 0-based; `POS == end == pam_positions`,
1149/1149 exact — no slack needed, ±2 bp over-joins).

**CAR-T gotchas:** join on **guide**, never `sample_name` — gold says `NS0011-ABTB1`, the pipeline
says `ABTB1-KO-DNA`, and joining on sample name silently drops 15 of 25 gold samples. Reuse
`guide_from_sample` + `GUIDE_ALIAS` from `validate_recall.py`. The `guide` column in
`wgs_hotspot_scores.csv` is a **recurrence proxy from `add_recurrence()`**, not a join key — pass
`--training` or `--samplesheet` or the denominator silently empties. Truth-side repairs needed:
sheet typo `CTLA41`→`CTLA4`, and `CART_NS0011-CREBRF.edited_reads.xlsm` has a corrupted `chrom`
cell (`"c"` for chr5:173,090,439).

**`is_target` is not truth** — 10 on-target sites were reviewed and rejected.

### Rebuilding the curated WGS label

`bin/build_curated_wgs_label.py`. The manual review was **not** positive-unlabeled: the process was
take every indel with `indel_fraction >= 0.05` **AND** `indel_reads >= 10`, then adjudicate every
row by eye. So **inside that stratum a blank `manual_review` means REJECTED**, and those rejects are
the best negatives available — they fooled a rules filter but not a human. Confirmed by the data:
among the 53 confirmed WGS edits the minimum `indel_fraction` is 0.0561 and the minimum
`indel_reads` is 11, both just inside the stated cut-offs. `indel_count >= 10` is the **wrong**
column (loses 4 confirmed) — `indel_count` is allelic diversity, `indel_reads` is read support.

### Numbers safe to quote

- **CAR-T recall 52/52 = 1.000** vs manual review (49/49 on-target, 3/3 off-target). Model score
  alone gives 47/52; the **high-evidence rescue supplies the last 6**.
- **Off-target review burden falls 20.8×** (83 gated → 4 kept) **with zero on-target loss**
  (77 gated → 77 kept). *All 79 drops are off-target rows.*
- **Precision 0.889** — but say "production configuration", and see the caveat below.

### Numbers that are traps

- **Do NOT quote "99,238 → 81".** That credits the five rules with the evidence gate's work. The
  gate (`bin/review_filter.py:334`, `indel_reads >= min_reads & indel_fraction >= min_vaf`) does 99,238 → 160 by
  itself; the rules only ever see the 160. That gate is *exactly the filter the lab already applied
  by hand in Excel*, so the honest baseline is the post-gate count — **160** for the 32-sample
  cohort, **238** for the 25-sample one. "238 reviewed by hand" is defensible; "238 candidates" is
  not.
- **0.877 vs 0.889 is not a conflict** — different rows of the `docs/NOISE_MODEL_EXPERIMENT.md` arm
  table. 0.877 = the PoN-only / no-cohort-equivalence arm; **0.889 = production**. Never quote
  either bare. Retained edits are 64/64 in every arm.
- **Precision 0.889 predates the current caller.** It was measured on `results_cart_ponfix`
  (2026-08-13) where the same 99,238 rows and the same gate passed **1,498** rows; the current
  caller passes **479**. `find_edited_reads.py` changed. **It has not been re-adjudicated against
  the current 86-row queue.** Say so whenever quoting it.
- **Depth is ~202–205×, not 60×.** The "16% VAF floor at 61×" belongs to a *different* dataset
  (the WGS negative control) and must not be quoted for the CAR-T cohort. At 202× the pipeline's
  2-read/0.5% gate is an effective floor of **0.99% VAF**. The lab's historical 10-read hand filter
  is **4.95%** — *that* is where "5% detector" comes from, and it describes the old manual
  threshold, not this pipeline.
- **A BND queue row is not an event.** The caller reports each junction from both ends, at ±4 bp
  jitter, under both orientations. The CAR-T cohort's **25 rows are 8 junctions**. Never quote a row
  count as an event count, and never quote a row's `reads` as event support — ARID4A's rows say 3–4
  each; the junction carries 18.
- **AAVS1's "8 real edits" are 2 distinct sites seen once per sample**, not 8 sites.

---

## 5. The breakend result (corrected 2026-08-20)

The PI uses "breakend" to mean *translocation*. The answer to "are any of these translocations?" is
**no — all 8 are intra-chromosomal.**

Checking that exposed a shipped labelling bug: `review_filter_bnd.py`'s classifier keyed on
`(is_target, far_end_on_target, interchromosomal)` and **never read `strands`**, so it called
everything "multi-cut deletion". Every queue row is `+-`/`-+` — the **inverted** adjacency.
**The 8 junctions are inversions.**

Two independent confirmations: (1) all 8 events carry **both junctions of the reciprocal pair**
(support balance 0.43–0.80), which a deletion cannot produce — it makes one junction; (2) across
every raw breakend record at the 8 cut pairs, unfiltered, there are **0 collinear junctions vs 5–16
inverted**, so the deletion product is absent, not sub-threshold. Structural reason the old label
was wrong by construction: a same-chromosome, **same-strand** split read is resolved to DEL/DUP/INS
at `find_edited_reads.py:420` and never becomes a BND — so a same-chromosome BND can *only* be
inverted.

**The coverage test is underpowered — do not cite it.** At 5–20% allele fraction the expected drop
is 3–10% and flank scatter is ±10%; it came back 0.87–1.01, uninformative.

**Biology:** 7 of 8 target files carry two different guides against the KO gene (BRAF has one), so
these are dual-guide events where the intervening fragment was **flipped and re-ligated** rather
than lost. **No junction anywhere touches TRAC or TRBC** despite both being present in all 8 target
files — that is a real negative result about the two guides common to every construct.

**Still open:** `bin/bnd_from_indels_to_vcf.py` was fixed (it wrote **zero records in every run
ever** — doubled escapes meant the header never split) and now emits **1,022 records across 32
samples**. But **that does not populate the HTML report, and I said at one point that it did.**
`bin/make_scge_report.qmd` has **zero references to `bnd_vcf`**; the report's SV panel renders
`tables.on_target_sv_transgene`, an unrelated transgene-junction table that is empty in all 32
samples. Two independent defects; only the VCF one is fixed. **Surfacing breakends in the report is
a QMD change with the data already waiting** — this is the single highest-value small task left.

---

## 6. The noise model / panels

Full treatment in `docs/NOISE_MODEL.md`, `docs/NOISE_MODEL_ASSUMPTIONS.md`,
`docs/NOISE_MODEL_VALIDATION.md`, `docs/PANEL_AS_FILTER.md`. Operational points only:

- **The matched-control beta-binomial reproduces the panel-of-normals exactly** with no cohort
  required. The PoN is gone.
- **`AQ_MIN` is 5, not 30.** The depth floor rescales AQ; AQ≥10 starts costing confirmed edits.
- **The clean-control trap:** a control with 0 alt reads at depth *d* bounds background at ~1/*d*,
  not 0. `apply_depth_floor` handles it. Currently masked by the `VAF>=0.005` gate — **it stops
  being masked if that gate is ever lowered.**
- **Omitting `--pon` does not disable rule 4** — it falls back to cross-guide recurrence. Use
  `--max-guides 9999` to actually turn it off.
- **`review_filter.py`'s own defaults are 10 reads / 5% VAF, but the pipeline passes 2 / 0.005.**
  Any offline rerun must pass them explicitly or it will not reproduce the baseline.
- **Cross-guide recurrence needs all guides in ONE invocation.** Single-guide runs silently get no
  rule 4. It now warns; `--strict-fallback` makes it fatal.

### Two different params, only one reaches rule 4 — this caused real confusion

- **`params.sv_noisefile`** (`nextflow.config:213`) → **DRAGEN** `--sv-systematic-noise`. The
  breakend arm never sees it.
- **`params.review_sv_noise`** (`nextflow.config:172`) → **rule 4**, default **null = OFF**.

⇒ `dropped, DRAGEN systematic noise : 0` in a run means **the rule did not run**, not that nothing
matched.

**`sv_noisefile` stays on IDPF v3.0.0. This is the PI's decision and the question is CLOSED.** It
is harmless for everything in these docs — the breakend caller reads SA-tag split reads only, never
DRAGEN SV calls. It only shapes DRAGEN's own `*.sv.annotated.vcf.gz` and the report SV table
(~56% of passing calls suppressed).

**The real risk to guard: do not wire IDPF into `review_sv_noise`.** It flags **25/25 real
junctions** at every slop and would erase the entire breakend finding. IDPF covers 62.7% of the
genome, so its 85.9% flag rate is exactly its 86.1% coverage null (1.00×) — and it points
*backwards*, flagging on-target real edits more than off-target (OR 0.45, CI excludes 1). Warning
sits at `nextflow.config:160-171`.

**The SNV/indel panel is the opposite case and does work.** `params.snv_noisefile`
(`IDPF_WGS_hg38_v.2.0.0_systematic_noise.snv.bed.gz`) is a *different file*. On AAVS1: 443
candidates clear the gate, **54 flagged (12.2%), all 54 off-target, 0 confirmed edits flagged**,
212 left to review. It works because it is **specific** — only 2.05M of 97.85M records qualify
(≥3 donors + an indel-capable allele) = **≤0.332% of the genome**, so 12.2% is ~37× the null.
One-liner: *a noise filter is only useful if it is specific — ours is on indels, the SV one is not.*

---

## 7. Open threads, ranked

1. **Open the PR.** Nothing else on this list matters if the branch is never merged.
2. **Surface breakends in the HTML report** (§5). Small QMD change, data already produced.
3. **Retire the pytest failures** by running the suite inside
   `ghcr.io/dhslab/docker-scge-offtarget:260710`. The missing module is `pyranges`.
4. **The model does not beat plain rules on AAVS1** (precision 0.067 vs 0.095). The 40× review-queue
   win is real but it comes from **`cut_dist`, which is a rule**, not from the model. Be honest
   about this.
5. **The ceiling is the label, not the features or the model.** Swapping only the training label
   (7 features held fixed) costs 6× precision. A 10-feature curated model improves AAVS1 precision
   2.1× (28→12 FPs) **but regresses CAR-T recall 52/52 → 47/52**, so it was *not* shipped —
   `assets/models/wgs_shape_model.pkl` is still v1. The 5 misses all sit below the rescue's
   `indel_frac >= 0.15` floor. Two concrete next experiments: (a) draw positives from the ECS review
   (all 56, including the 7 below WGS 0.05) with negatives from the WGS stratum — the current label
   draws both classes from the WGS stratum and under-represents low-VAF positives; (b) lower or
   depth-adapt the rescue floor. Also try fewer features given only 69 positives.
   **A CNN is the wrong move** — more capacity fits the bad label better and performs worse.
6. **The hard false-positive set.** GM24385-unedited scored at AAVS1 hotspots gives **175 certain
   FPs** (3.44%), no review needed. Of the 21 low-VAF ones, **7 survive every rule we have** —
   chr6:29918619, chr7:73852004, chr8:10469359, chr8:142037418, chr11:113381806, chr12:116647121,
   chr22:48757172. That is the model's actual job and its evaluation set. (chr22:48757172 has 18
   distinct indel lengths at cut_dist 8 with 2,769 reads — indistinguishable from a real edit by any
   rule.) Note this also **tempers the `cut_dist` claim**: against a proper negative control the
   median is 14 bp and **40.6% are within 10 bp**, so `cut_dist<=10` removes 59%, not the ~82% seen
   against the curated TN set.
7. **The well-posed WGS sensitivity target, label-free:** reaching real 5% sensitivity in WGS needs
   `reads>=3`, which costs **9 → 34 artifact FPs/sample**. The model's job is to win that ~4× back
   at fixed sensitivity. Measurable with no manual review.
8. **The uncorrected-ECS rerun is verified GO but unexecuted.** See §8.

---

## 8. The uncorrected-ECS rerun (verified, not run)

Realign the AAVS1 ECS capture with UMI consensus disabled, to get a noisy WGS-like training set
labelled from the corrected run. **Verified 2026-07-31: the SLURM stub passes and the diff is
exactly one flag.**

- Route: **AWS Batch via `dragenflow`**, `-profile dragenaws,alignonly,idtumi -c slurm.config
  --readfamilysize 1`. Canonical checkout
  `/storage2/fs1/dspencer/Active/spencerlab/abonney/git_runs/dragenflow`.
- Queue `dragen-queue_v4-4-6` **tested end-to-end from the SLURM node** — job SUCCEEDED, container
  reports `dragen Version 13.021.779.4.4.6`, byte-identical to the corrected run. Only `-profile ris`
  (LSF; `bsub` does not exist here) must become SLURM.
- **Pass `cpus 24` / `memory 240.GB` explicitly** or jobs stick in RUNNABLE forever
  (`MISCONFIGURATION:JOB_RESOURCE_REQUIREMENT` — the job definition's default 256000 MB does not fit
  an f2.6xlarge once ECS overhead is taken).
- **Use `readfamilysize = 1`, NOT `umi = null`.** `params.umi` is overloaded: when it is set the
  pipeline runs fastp for adapter trimming; when null, DRAGEN does it instead. Unsetting `umi`
  silently swaps trimmers *and* leaves ~3× PCR duplicates uncollapsed. `readfamilysize=1` changes
  exactly one flag: `--umi-min-supporting-reads 3` → `1`.
- **Cost is not the risk** (~$12 incremental for the original Jan run; ceiling ~$30–50). **Staging
  is**: a real run uploads ~853 GB to S3, and fastp writes a similar volume locally first (~900 GB
  scratch). Do not test the queue with a full `-stub-run` — submit a trivial job directly.
- **Data gain is ~3×, not 15×.** 7.4× in raw reads, but independent molecules (families) are
  300,363,223 vs 99,904,287 emitted = **3.0×**, and only the molecule ratio bounds independent draws.

---

## 9. Assets and reference data

- Container: `ghcr.io/dhslab/docker-scge-offtarget:260710` (matplotlib 3.10.8, pysam 0.24.0,
  sklearn 1.8.0 — **the pin must match the pickle**). `docker-scge:latest` has the caller deps but
  **no matplotlib**, which is why figures run in the other image.
- `ghcr.io/dhslab/docker-casoffinder-bulge:latest` — pushed and working.
- Lab images all live in **one repo, `dhslab/dhslab-docker-images`**, one per `docker-*/` dir;
  a push to main auto-builds any changed `**/Dockerfile` to `ghcr.io/dhslab/<dirname>` with
  `:latest` + `:YYMMDD`. **The directory name becomes the image name** the modules pin.
- CRISPRme prebuilt index (6.8 G, ~17 min to build):
  `/storage2/fs1/dspencer/Active/clinseq/projects/scge/data/refdata/crisprme_hg38/`. **Reuse gotcha:** `complete-search` resolves `genome_library`
  relative to **CWD**, not to `--genome` — an absolute `--genome` from a foreign CWD silently
  *rebuilds* the index. `crisprme.nf` symlinks `Genome` and `genome_library` into the task CWD.
- Cas-OFFinder **3.0.0 does bulges natively**; the 2016 `cas-offinder-bulge` wrapper is obsolete and
  crashes on v3. Call `cas-offinder <input> C <out>` directly. Also: Nextflow launches apptainer with
  `--no-home`, so pocl cannot write its kernel cache and dies with "No OpenCL devices found" — fixed
  by `export POCL_CACHE_DIR="${PWD}/.pocl_cache"` in the module.
- DRAGEN noise files: `/storage2/fs1/dspencer/Active/spencerlab/refdata/hg38/dragenfiles` (this is
  the one `nextflow.config` points at; a second copy exists under `Active/shared/refdata/hg38/`).
- References differ by cohort: AAVS1 ECS used `hg38_mgi_patch.fa` + refdir
  `dragen_hg38_cg_rna_cnv_v4.4.6`; the CAR-T/WGS side uses `hg38_PLVM_CD19_CARv4_cd34.fa`.
  **Do not mix them.**

---

## 10. Things I got wrong, so you don't re-derive them

- I said the 8 junctions were deletions. **They are inversions** (§5).
- I said the SV VCF fix populated the HTML report. **It does not** (§5).
- I said unedited WGS controls did not exist. **They do** — DRAGEN tumor/normal writes the normal to
  `<prefix>.cram` and the tumor to `<prefix>_tumor.cram`, named after the output prefix rather than
  its contents. Three donors: `CD34-CART-DNA` (21 samples), `CART_NS0027-unedited` (8),
  `CART_NS0065-unedited` (3).
- I concluded "assumption 8 fails" in the noise-model audit. **That was a selection artifact** from
  pooling 3% of rows selected for carrying a called indel. Assumption 8 holds where testable.
- The 11.7× clean-locus gap is **5 loci carrying donor-private germline** the caller had already
  suppressed — not a model defect.
- `control_indel_reads` is **not** constant-0, but it is untrustworthy: `--max-in-control` defaults
  to 5 (not 0), and the column is a floor-rounded *mean* over only variants that passed that filter,
  so it silently drops the strongest germline evidence. The per-event `control_alt_counts` (field 11
  of `indel_info`) **is** 0 by construction and carries no information.
- **The test that settles coordinate questions here:** run the caller with `--edited-bam` and
  `--control-bam` set to the **same file**. Every called event must then be control-supported, so
  anything escaping `-x` is a provable miss. Label-free, ~3 min on a 92-site subset. Two careful
  code readings gave the wrong answer on this; this test gave the right one first try.
