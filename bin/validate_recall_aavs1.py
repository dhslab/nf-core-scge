#!/usr/bin/env python3
"""
validate_recall_aavs1.py — WGS performance against the CURATED AAVS1 truth set.

Why this is separate from validate_recall.py
--------------------------------------------
`validate_recall.py` scores the CART cohort: its `--gold` is a merged review workbook keyed by
`sample_name`/`manual_review`, joined on the guide, and it is deliberately RECALL-ONLY because
that table has no confirmed negatives (NaN there means "not reviewed", not "rejected").

AAVS1's curation is a different artifact entirely — two read-level tables, one of confirmed
edit-supporting reads and one of confirmed non-edit reads:

    crispr_ml/AAVS1_training_tp.tsv    77,317 reads
    crispr_ml/AAVS1_training_tn.tsv    72,284 reads

Because it has a real negative class, precision and specificity ARE defined here, unlike CART.
Different schema, different join, different metrics — hence a separate script rather than more
flags on the CART one.

What the truth set actually contains
------------------------------------
Collapsed to sites: 2 TP sites and 6,037 TN sites. Both TP sites (chr19:55115732 and
chr19:55115752) are **ON-TARGET** — AAVS1 site5 and site14. There are no human-confirmed
off-target edits in this cohort, so "100% recall" here means both on-target edits were recovered
from WGS alone. Say that out loud rather than letting it be discovered.

Three traps, each of which silently produces a WRONG answer
-----------------------------------------------------------
1. The TSVs are RAGGED: the header names 4 columns but each row carries a full SAM record across
   many unnamed ones. `usecols=['Label','VCF_Chrom','VCF_Pos']` returns an EMPTY frame — which
   reads downstream as "0 TP sites", i.e. a passing-looking 0/0, not an error. Use positional
   `usecols=[0,1,2]`.
2. These are READ-level labels. An edited site contains both edit-supporting reads AND ordinary
   reference reads, so **both TP sites also appear in the TN table**. Subtract the TP sites from
   the TN set or they score as 2 phantom false positives.
3. `VCF_Pos` is 1-based; the pipeline's `start` is a 0-based BED-style coordinate. Join with a
   small slack (default +/-2) or every site misses.

Usage
-----
  validate_recall_aavs1.py \
      --training results_offtarget_aavs1_igv/offtarget/training.tsv \
      --tp /storage2/.../crispr_ml/AAVS1_training_tp.tsv \
      --tn /storage2/.../crispr_ml/AAVS1_training_tn.tsv \
      --out aavs1_audit.csv --require-recall 1.0
"""
import argparse
import sys

import pandas as pd

DEFAULT_TP = ("/storage2/fs1/dspencer/Active/clinseq/projects/scge/"
              "crispr_ml/AAVS1_training_tp.tsv")
DEFAULT_TN = ("/storage2/fs1/dspencer/Active/clinseq/projects/scge/"
              "crispr_ml/AAVS1_training_tn.tsv")
DETECT_PATTERN = "LIKELY EDIT"
SLACK = 2


def read_curated(path, what):
    """Read a read-level curation TSV -> distinct (chrom, pos) sites.

    Positional usecols is deliberate: see trap 1 in the module docstring.
    """
    d = pd.read_csv(path, sep="\t", usecols=[0, 1, 2], header=0,
                    names=["label", "chrom", "pos"], low_memory=False)
    n_raw = len(d)
    # the ragged tail can smear stray values into these columns; keep only real loci
    d = d[d["chrom"].astype(str).str.match(r"^chr")]
    d["pos"] = pd.to_numeric(d["pos"], errors="coerce")
    d = d.dropna(subset=["pos"])
    d["pos"] = d["pos"].astype(int)
    sites = d.drop_duplicates(["chrom", "pos"])[["chrom", "pos"]].reset_index(drop=True)
    if sites.empty:
        sys.exit(f"ERROR: {what} ({path}) yielded 0 sites from {n_raw} rows. The file is ragged; "
                 "this is what named usecols does. Read columns positionally.")
    print(f"  {what:<3} {n_raw:>7,d} reads -> {len(sites):>5,d} distinct sites")
    return sites


def load_truth(tp_path, tn_path):
    """Curated site sets, with the TP sites removed from the TN set (trap 2)."""
    print("curated AAVS1 truth:")
    tp = read_curated(tp_path, "TP")
    tn = read_curated(tn_path, "TN")
    tp_set = {(c, p) for c, p in zip(tp["chrom"], tp["pos"])}
    before = len(tn)
    tn = tn[[(c, p) not in tp_set for c, p in zip(tn["chrom"], tn["pos"])]]
    dropped = before - len(tn)
    if dropped:
        print(f"  removed {dropped} edit-site(s) from the TN set (a real edit also carries "
              f"reference reads, so it appears in both tables)")
    return tp_set, {(c, p) for c, p in zip(tn["chrom"], tn["pos"])}


def _near(chrom, pos, site_set, slack=SLACK):
    """1-based VCF vs 0-based pipeline start (trap 3)."""
    return any((chrom, pos + off) in site_set for off in range(-slack, slack + 1))


def label_frame(df, tp_set, tn_set, slack=SLACK):
    """Annotate a scored table with is_TP / is_TN / detected. Reused by compare_to_baseline."""
    out = df.copy()
    out["detected"] = out["verdict"].astype(str).str.contains(DETECT_PATTERN, na=False)
    out["is_TP"] = [_near(c, p, tp_set, slack) for c, p in zip(out["chrom"], out["start"])]
    out["is_TN"] = [_near(c, p, tn_set, slack) for c, p in zip(out["chrom"], out["start"])]
    return out


def confusion(df):
    tp = int(df.loc[df["is_TP"], "detected"].sum())
    fn = int((~df.loc[df["is_TP"], "detected"]).sum())
    fp = int(df.loc[df["is_TN"], "detected"].sum())
    tn = int((~df.loc[df["is_TN"], "detected"]).sum())
    return tp, fn, fp, tn


def report(tp, fn, fp, tn):
    recall = tp / (tp + fn) if (tp + fn) else float("nan")
    spec = tn / (tn + fp) if (tn + fp) else float("nan")
    prec = tp / (tp + fp) if (tp + fp) else float("nan")
    print()
    print("=" * 70)
    print("WGS-ONLY PERFORMANCE vs CURATED AAVS1 TRUTH")
    print("=" * 70)
    print(f"  TP {tp}   FN {fn}   FP {fp}   TN {tn}")
    print(f"  recall      = {tp}/{tp + fn} = {recall:.3f}")
    print(f"  specificity = {tn}/{tn + fp} = {spec:.4f}")
    print(f"  precision   = {tp}/{tp + fp} = {prec:.3f}")
    print()
    print("  NOTE: both curated TP sites are ON-TARGET (AAVS1 site5 / site14). This cohort has")
    print("        no human-confirmed off-target edits, so recall here is recovery of the")
    print("        on-target edits from WGS alone.")
    print("  NOTE: precision is against CURATED negatives. It is NOT comparable to the")
    print("        precision in offtarget_metrics.txt, which uses the ECS label.")
    return recall


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--training", required=True, help="training.tsv from the OFFTARGET run")
    ap.add_argument("--tp", default=DEFAULT_TP, help="AAVS1_training_tp.tsv")
    ap.add_argument("--tn", default=DEFAULT_TN, help="AAVS1_training_tn.tsv")
    ap.add_argument("--slack", type=int, default=SLACK,
                    help="bp slack for the 1-based/0-based join (default 2)")
    ap.add_argument("--max-cut-dist", type=float, default=None,
                    help="also report the confusion matrix with a |modal_pos - start| gate "
                         "applied, i.e. require the observed indel to sit near the predicted cut")
    ap.add_argument("--out", help="write the per-site audit table here")
    ap.add_argument("--require-recall", type=float,
                    help="exit non-zero if recall falls below this (for CI)")
    args = ap.parse_args()

    tp_set, tn_set = load_truth(args.tp, args.tn)
    df = pd.read_csv(args.training, sep="\t", low_memory=False)
    lab = label_frame(df, tp_set, tn_set, args.slack)

    matched_tp = int(lab["is_TP"].sum())
    if matched_tp == 0:
        sys.exit("ERROR: no curated TP site matched the scored table. The join is broken — check "
                 "the chrom naming and --slack. An empty denominator is not a pass.")
    print(f"  matched in scored table: {matched_tp} TP row(s), {int(lab['is_TN'].sum())} TN row(s)")

    tp, fn, fp, tn = confusion(lab)
    recall = report(tp, fn, fp, tn)

    print()
    print("--- the curated TP sites ---")
    cols = [c for c in ("sample", "chrom", "start", "modal_pos", "indel_frac", "spanning",
                        "ctrl_if", "score", "verdict", "ecs_if") if c in lab.columns]
    print(lab.loc[lab["is_TP"], cols].to_string(index=False))

    if args.max_cut_dist is not None and {"modal_pos", "start"} <= set(lab.columns):
        g = lab.copy()
        g["cut_dist"] = (g["modal_pos"] - g["start"]).abs()
        # a site fails the gate if it was called but the indel sits far from the predicted cut
        g["detected"] = g["detected"] & (g["cut_dist"] <= args.max_cut_dist)
        print()
        print(f"### with cut-site gate |modal_pos - start| <= {args.max_cut_dist:g} ###")
        report(*confusion(g))

    if args.out:
        keep = lab[lab["is_TP"] | lab["is_TN"]].copy()
        if {"modal_pos", "start"} <= set(keep.columns):
            keep["cut_dist"] = (keep["modal_pos"] - keep["start"]).abs()
        keep["curated"] = ["TP" if t else "TN" for t in keep["is_TP"]]
        keep.to_csv(args.out, index=False)
        print(f"\nwrote per-site audit -> {args.out}  ({len(keep)} rows)")

    if args.require_recall is not None and recall < args.require_recall:
        sys.exit(f"FAIL: recall {recall:.3f} < required {args.require_recall:.3f}")


if __name__ == "__main__":
    main()
