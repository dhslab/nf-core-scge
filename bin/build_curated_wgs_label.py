#!/usr/bin/env python3
"""
build_curated_wgs_label.py — the training label the shape model should actually learn.

Why this exists
---------------
`join_training_table.py` derives `label` from the ECS indel threshold. On the CART panel that
yields 73,770 "positives" out of 99,308 rows, most of them ECS assay noise or germline. A model
fit to that target scores *worse* in reality even while its AUC against the target goes up
(measured: swapping only the label cost 6x precision on the AAVS1 curated truth).

The defensible label is the human one, and it is a genuine TWO-CLASS label -- which the code
previously assumed it was not. The manual WGS review worked like this:

  1. take every called indel with `indel_fraction >= 0.05` AND `indel_reads >= 10`
  2. adjudicate every row in that set by eye in IGV

So inside that stratum, a blank `manual_review` means **rejected**, not unreviewed. Those
rejects are the most valuable negatives available: sites that passed a rules filter and a human
still said no -- exactly the false positives the model exists to remove.

The thresholds are confirmed by the data rather than assumed: among the 53 confirmed WGS edits
the minimum `indel_fraction` is 0.0561 and the minimum `indel_reads` is 11, both just inside the
stated cut-offs, and no confirmed edit falls outside the stratum.

  gold_wgs, indel_fraction >= 0.05 & indel_reads >= 10
      -> 241 rows: 53 confirmed edits, 188 human-rejected negatives

Rows OUTSIDE the stratum are genuinely unreviewed and are dropped, not labelled 0.

Note the ECS review used a lower VAF floor (7 of its confirmed edits sit below 0.05), which is
reasonable at 5000x. Do not reuse these thresholds for the ECS gold.

Usage
-----
  build_curated_wgs_label.py --training cart_training_newfeats.tsv \
      --gold-xlsx ".../cart_wgs/cart_wgs_merged.xlsx" --out cart_training_curated.tsv
"""
import argparse
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
# The gold sheets key samples by review-sheet name (`NS0011-ABTB1`) while the pipeline keys
# them by CRAM sample (`ABTB1-KO-DNA`). Joining on sample_name silently drops 15 of 25 gold
# samples. Both carry the GUIDE, so the guide is the join key -- the same resolution
# validate_recall.py already implements, alias table included.
from validate_recall import guide_from_sample, GUIDE_ALIAS  # noqa: E402

CONFIRMED = {"1", "1.0"}
MIN_IF = 0.05
MIN_INDEL_READS = 10


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--training", required=True, help="training table carrying the WGS features")
    ap.add_argument("--gold-xlsx", required=True, help="cart_wgs_merged.xlsx")
    ap.add_argument("--sheet", default="gold_wgs")
    ap.add_argument("--min-if", type=float, default=MIN_IF)
    ap.add_argument("--min-indel-reads", type=int, default=MIN_INDEL_READS)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    g = pd.read_excel(args.gold_xlsx, sheet_name=args.sheet)
    g["chrom"] = g["chrom"].astype(str)
    g["start"] = pd.to_numeric(g["start"], errors="coerce").astype("Int64")
    conf = g["manual_review"].astype(str).str.strip().isin(CONFIRMED)
    iv = pd.to_numeric(g["indel_fraction"], errors="coerce")
    ir = pd.to_numeric(g["indel_reads"], errors="coerce")
    in_stratum = (iv >= args.min_if) & (ir >= args.min_indel_reads)

    # sanity: the stratum must contain every confirmed edit, or the thresholds are wrong
    lost = int((conf & ~in_stratum).sum())
    if lost:
        sys.exit(f"ERROR: {lost} confirmed edit(s) fall OUTSIDE the stratum "
                 f"(if>={args.min_if}, indel_reads>={args.min_indel_reads}). The review filter "
                 f"must be wrong -- a confirmed edit cannot be outside the reviewed set.")

    rev = g[in_stratum].copy()
    rev["curated_label"] = conf[in_stratum].astype(int)
    rev["guide"] = rev["sample_name"].map(lambda x: guide_from_sample(x, GUIDE_ALIAS))
    # one site can be reviewed in several replicates of the same guide; if any replicate was
    # confirmed the site is a real edit, so take the max
    keys = (rev.groupby(["guide", "chrom", "start"], as_index=False)["curated_label"].max())
    print(f"reviewed stratum: {len(keys)} sites over {keys.guide.nunique()} guides "
          f"({int(keys.curated_label.sum())} confirmed, "
          f"{int((keys.curated_label == 0).sum())} human-rejected)")

    t = pd.read_csv(args.training, sep="\t", low_memory=False)
    t["chrom"] = t["chrom"].astype(str)
    t["start"] = pd.to_numeric(t["start"], errors="coerce").astype("Int64")
    if "guide" not in t.columns:
        sys.exit("ERROR: training table has no `guide` column to join on.")
    m = t.merge(keys, on=["guide", "chrom", "start"], how="inner")
    if m.empty:
        sys.exit("ERROR: the gold/training join produced 0 rows -- check guide resolution.")
    matched = set(keys.guide) & set(t.guide)
    missing = sorted(set(keys.guide) - set(t.guide))
    print(f"  guides matched: {len(matched)}/{keys.guide.nunique()}"
          + (f"; absent from the scored cohort: {missing}" if missing else ""))

    # the curated label REPLACES the ECS-threshold label; keep the old one for comparison
    m = m.rename(columns={"label": "label_ecs_threshold"})
    m["label"] = m["curated_label"]
    m.to_csv(args.out, sep="\t", index=False)
    n_pos = int((m.label == 1).sum())
    print(f"wrote {args.out}: {len(m)} rows joined to features "
          f"({n_pos} positive / {len(m) - n_pos} negative) over {m['sample'].nunique()} samples")
    if "label_ecs_threshold" in m.columns:
        agree = int((m.label == m.label_ecs_threshold).sum())
        print(f"  the ECS-threshold label agrees with the human on {agree}/{len(m)} of these rows")


if __name__ == "__main__":
    main()
