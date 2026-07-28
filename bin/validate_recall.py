#!/usr/bin/env python3
"""
validate_recall.py — WGS recall against the HUMAN-REVIEWED gold standard.

Why this exists
---------------
`recall_vs_vaf.py` measures recall against the ECS label (`ecs_is_edit`), whose
denominator is set by `offtarget_ecs_edit_threshold`. At the default of 0.0 ANY nonzero
ECS indel fraction counts as an edit, so the denominator fills with ECS noise: in the
CART run 47,170 of ~64,000 "positives" sit below 0.5% VAF. Recall against that
denominator (~2%) says nothing about whether the pipeline finds real CRISPR edits.

The defensible denominator is the manual review: the per-sample `*.edited_reads.xlsm`
curation merged into `cart_ecs_merged.csv.gz`, where `manual_review == 1` marks a
human-confirmed real CRISPR edit. This script scores the pipeline against exactly that
set and prints the misses, so "100% recall" is a number anyone can re-derive rather than
a claim.

Usage
-----
  validate_recall.py --scores  results/offtarget/wgs_hotspot_scores.csv \
                     --gold   "Manual Indel Review/cart_ecs/cart_ecs_merged.csv.gz"

Joins on (guide, chrom, start). The gold table keys samples by their review-sheet name
(`CART_NS0011-ABTB1`); the pipeline keys them by CRAM sample (`ABTB1-KO-DNA`). Both
carry the guide, so the guide is the join key and `--guide-alias` patches the handful of
sheet-name typos that would otherwise silently drop truth rows.
"""
import argparse
import re
import sys

import pandas as pd

# Sheet-name typos in the manual-review workbooks that do not match the pipeline's guide
# names. These are DATA fixes, not logic: without them the affected truth rows join to
# nothing and inflate recall by vanishing from the denominator.
GUIDE_ALIAS = {"CTLA41": "CTLA4"}

# A single corrupted cell: CART_NS0011-CREBRF.edited_reads.xlsm carries chrom "c" for
# chr5:173,090,439 (CREBRF is chr5q35.1). Same reasoning — repair, do not drop.
CHROM_REPAIR = {"c": "chr5"}

CONFIRMED = {"1", "1.0", "1?"}


def guide_from_sample(name, alias):
    """`CART_NS0011-ABTB1` -> ABTB1; `CART_NS0027-CTLA41_2` -> CTLA4."""
    g = str(name).split("-", 1)[1] if "-" in str(name) else str(name)
    g = re.sub(r"_\d+$", "", g)
    return alias.get(g, g)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scores", required=True, help="wgs_hotspot_scores.csv from score.py")
    ap.add_argument("--gold", required=True,
                    help="merged manual-review table with a manual_review column")
    ap.add_argument("--samplesheet",
                    help="offtarget samplesheet (sample,guide,...) — the AUTHORITATIVE "
                         "sample->guide map. Without it the scored table's own `guide` "
                         "column is used, which add_recurrence() derives from the sample "
                         "name and is only a recurrence proxy (e.g. 'ABTB1-KO-DNA').")
    ap.add_argument("--training",
                    help="training.tsv — alternative sample->guide source (it carries the "
                         "samplesheet-resolved guide)")
    ap.add_argument("--guide-alias", action="append", default=[], metavar="SHEET=PIPELINE",
                    help="extra guide-name alias; repeatable")
    ap.add_argument("--hi", type=float, default=0.60,
                    help="model-score threshold, for the score-only comparison")
    ap.add_argument("--out", help="write the per-site audit table here")
    ap.add_argument("--require-recall", type=float,
                    help="exit non-zero if recall falls below this (for CI)")
    args = ap.parse_args()

    alias = dict(GUIDE_ALIAS)
    for a in args.guide_alias:
        k, _, v = a.partition("=")
        alias[k] = v

    gold = pd.read_csv(args.gold, low_memory=False)
    if "manual_review" not in gold.columns:
        sys.exit("ERROR: --gold has no manual_review column; this is not a review table")
    gold = gold[gold["manual_review"].astype(str).str.strip().isin(CONFIRMED)].copy()
    gold["guide"] = gold["sample_name"].map(lambda s: guide_from_sample(s, alias))
    gold["chrom"] = gold["chrom"].astype(str).replace(CHROM_REPAIR)
    gold["start"] = pd.to_numeric(gold["start"], errors="coerce").astype("Int64")
    gold = gold.rename(columns={"indel_fraction": "ecs_if"})

    sc = pd.read_csv(args.scores, low_memory=False)
    # Resolve sample -> guide. The scored table's own `guide` column comes from
    # add_recurrence(), which just strips a _1/_2 replicate suffix off the sample name —
    # fine for the recurrence grouping it was built for, wrong as a join key against the
    # review sheets. Prefer a real map when one is supplied.
    smap = None
    if args.samplesheet:
        ss = pd.read_csv(args.samplesheet)
        ss.columns = [c.strip().lower() for c in ss.columns]
        smap = dict(zip(ss["sample"], ss["guide"]))
    elif args.training:
        tt = pd.read_csv(args.training, sep="\t")
        smap = dict(zip(tt["sample"], tt["guide"]))
    if smap:
        sc["guide"] = sc["sample"].map(smap)
    elif "guide" not in sc.columns:
        sys.exit("ERROR: --scores has no guide column; pass --samplesheet or --training")
    sc["chrom"] = sc["chrom"].astype(str)
    sc["start"] = pd.to_numeric(sc["start"], errors="coerce").astype("Int64")
    sc["score"] = pd.to_numeric(sc["score"], errors="coerce")
    sc["_called"] = sc["verdict"].astype(str).str.contains("LIKELY EDIT", na=False)

    scored_guides = set(sc["guide"].dropna())
    unmatched = sorted(set(gold["guide"]) - scored_guides)
    if unmatched:
        print(f"NOTE: {len(unmatched)} gold guide(s) absent from the scored cohort "
              f"(not run / no CRAM), excluded from the denominator: {unmatched}")
        gold = gold[gold["guide"].isin(scored_guides)]
    # An empty denominator means the join key is broken, NOT that the pipeline passed.
    # Fail loudly — a silent 0/0 here would read as success.
    if gold.empty:
        sys.exit("ERROR: no gold guide matched the scored cohort — the sample->guide join "
                 "is broken. Pass --samplesheet (or --training) so guides resolve, and "
                 "check --guide-alias for review-sheet name typos.")

    if "call_basis" not in sc.columns:      # pre-rescue scored table
        sc["call_basis"] = ""
    m = gold.merge(sc, on=["guide", "chrom", "start"], how="left", suffixes=("_gold", ""))
    site = (m.groupby(["guide", "chrom", "start"], dropna=False)
              .agg(called=("_called", "max"), best=("score", "max"),
                   evaluated=("score", "count"),
                   verdict=("verdict", lambda s: "|".join(sorted(set(s.dropna().astype(str))))),
                   basis=("call_basis", lambda s: "|".join(sorted(set(s.dropna().astype(str)) - {""}))),
                   ecs_if=("ecs_if", "first"), is_target=("is_target_gold", "first"),
                   indel_frac=("indel_frac", "max"), ctrl_if=("ctrl_if", "max"),
                   conc_ratio=("conc_ratio", "max"), spanning=("spanning", "max"))
              .reset_index())
    site["called"] = site["called"].fillna(False).astype(bool)

    n = len(site)
    n_called = int(site["called"].sum())
    n_scoreonly = int((site["best"] >= args.hi).sum())
    n_rescued = int(site["basis"].astype(str).str.contains("high-evidence").sum())
    n_noteval = int((site["evaluated"] == 0).sum())

    print(f"\n== WGS recall vs manual review ({args.gold.split('/')[-1]}) ==")
    print(f"  human-confirmed edits in scored guides : {n}")
    print(f"  evaluated by the WGS arm               : {n - n_noteval}/{n}")
    print(f"  RECALL (verdict = LIKELY EDIT)         : {n_called}/{n} = {n_called / n:.3f}")
    print(f"  recall on model score >= {args.hi} alone     : {n_scoreonly}/{n} = "
          f"{n_scoreonly / n:.3f}")
    print(f"  recovered by the high-evidence rescue  : {n_rescued}")
    on = site[site["is_target"] == 1]
    off = site[site["is_target"] != 1]
    if len(on):
        print(f"  on-target  : {int(on['called'].sum())}/{len(on)}")
    if len(off):
        print(f"  off-target : {int(off['called'].sum())}/{len(off)}")

    miss = site[~site["called"]]
    if len(miss):
        print(f"\n-- {len(miss)} MISSED confirmed edit(s) --")
        with pd.option_context("display.width", 220, "display.max_rows", None):
            print(miss.drop(columns=["called"]).to_string(index=False))
    else:
        print("\n-- no missed confirmed edits --")

    if args.out:
        site.to_csv(args.out, index=False)
        print(f"\nwrote {args.out} ({n} sites)")

    if args.require_recall is not None and (n_called / n) < args.require_recall:
        sys.exit(f"FAIL: recall {n_called / n:.3f} < required {args.require_recall}")


if __name__ == "__main__":
    main()
