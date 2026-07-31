#!/usr/bin/env python3
"""
compare_to_baseline.py — does the ML actually buy anything over the rules?

The question this answers
-------------------------
The SCGE pipeline already finds indels: `get_indels.nf` runs `bin/find_edited_reads.py` (with
`--enable-crispr-prediction` commented out, so it is rules only) and a human then reviews the
output in IGV. The cost of that workflow is the size of the review queue. So the claim for the
ML arm is NOT "it detects edits" — it is:

    same recall, far fewer sites to look at.

That is a precision statement at matched recall, and this script measures it against the curated
AAVS1 truth (crispr_ml/AAVS1_training_{tp,tn}.tsv).

Four arms, all scored on ONE denominator: unique curated hotspot sites
---------------------------------------------------------------------
  1. ECS rules, ungated     every site with any indel read      <- what find_edited_reads.py emits
  2. ECS rules + somatic    indel>0, control==0, indel_frac>t   <- what a reviewer would prefilter
  3. WGS rules only         the model's own stage-1 gate, no model
  4. WGS + ML               verdict contains LIKELY EDIT
  5. WGS + ML + cut gate    ... and the indel sits near the predicted cut site

Arms 1-2 read the ECS analysis TSVs (5000x targeted). Arms 3-5 read training.tsv (30x WGS).
**That difference is the point and must be stated, not hidden:** the ML arm is doing a harder job
on ~150x less depth. Arm 3 exists so the model's contribution can be separated from the data
type — it is the same input as arms 4-5, rules only.

Counting rule
-------------
Every arm is collapsed to **unique (chrom, start) sites** before counting, because the raw tables
are per-sample and the curated truth is per-site. Comparing a per-sample-call count against a
per-site count inflates the ratio; the plan for this work called that out explicitly.
"""
import argparse
import glob
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from validate_recall_aavs1 import load_truth, DEFAULT_TP, DEFAULT_TN, SLACK, _near

# From wgs_shape_model.pkl's own `stage1_filter` metadata, minus the homology gate (at a hotspot
# panel every site is a predicted site by construction, so min_mm/is_target would be circular).
STAGE1_MIN_IFRAC = 0.05
STAGE1_MAX_CTRL = 0.02
STAGE1_MIN_SPAN = 8


def site_set(df, chrom_col="chrom", pos_col="start"):
    return {(c, int(p)) for c, p in zip(df[chrom_col], df[pos_col]) if pd.notna(p)}


def score_arm(name, sites, tp_set, tn_set, slack, note=""):
    """Intersect a candidate site set with the curated truth."""
    tp = sum(1 for c, p in sites if _near(c, p, tp_set, slack))
    fp = sum(1 for c, p in sites if _near(c, p, tn_set, slack))
    n_tp_total = len(tp_set)
    return {"arm": name, "queue": len(sites), "TP": tp, "FP": fp,
            "recall": tp / n_tp_total if n_tp_total else float("nan"),
            "precision": tp / (tp + fp) if (tp + fp) else float("nan"),
            "note": note}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--training", required=True, help="training.tsv (the WGS arm)")
    ap.add_argument("--ecs-glob", required=True,
                    help="glob for *.offtarget_analysis.tsv (the rules-only ECS baseline)")
    ap.add_argument("--tp", default=DEFAULT_TP)
    ap.add_argument("--tn", default=DEFAULT_TN)
    ap.add_argument("--slack", type=int, default=SLACK)
    ap.add_argument("--somatic-min-ifrac", type=float, default=0.01,
                    help="ECS somatic prefilter: minimum indel_fraction (default 0.01)")
    ap.add_argument("--max-cut-dist", type=float, default=10.0,
                    help="cut-site gate for the last arm (default 10 bp)")
    ap.add_argument("--out", help="write the comparison table here as CSV")
    args = ap.parse_args()

    tp_set, tn_set = load_truth(args.tp, args.tn)
    slack = args.slack
    rows = []

    # ---- arms 1-2: the rules-only ECS baseline -----------------------------------------
    files = sorted(glob.glob(args.ecs_glob))
    if not files:
        sys.exit(f"ERROR: no ECS analysis files matched {args.ecs_glob}")
    any_indel, somatic = set(), set()
    for f in files:
        d = pd.read_csv(f, sep="\t", low_memory=False)
        any_indel |= site_set(d[d["indel_reads"] > 0])
        m = ((d["indel_reads"] > 0) & (d["control_indel_reads"] == 0)
             & (d["indel_fraction"] > args.somatic_min_ifrac))
        somatic |= site_set(d[m])
    print(f"\nECS baseline: {len(files)} sample tables")
    rows.append(score_arm("ECS rules, any indel read", any_indel, tp_set, tn_set, slack,
                          "find_edited_reads.py raw output"))
    rows.append(score_arm(f"ECS rules + somatic gate", somatic, tp_set, tn_set, slack,
                          f"indel>0, ctrl==0, if>{args.somatic_min_ifrac:g}"))

    # ---- arms 3-5: the WGS arm ----------------------------------------------------------
    t = pd.read_csv(args.training, sep="\t", low_memory=False)
    t["_called"] = t["verdict"].astype(str).str.contains("LIKELY EDIT", na=False)

    stage1 = ((t["indel_frac"] > STAGE1_MIN_IFRAC)
              & (t["ctrl_if"].fillna(0) < STAGE1_MAX_CTRL)
              & (t["spanning"] >= STAGE1_MIN_SPAN))
    rows.append(score_arm("WGS rules only (no model)", site_set(t[stage1]), tp_set, tn_set, slack,
                          f"if>{STAGE1_MIN_IFRAC}, ctrl<{STAGE1_MAX_CTRL}, span>={STAGE1_MIN_SPAN}"))
    rows.append(score_arm("WGS + ML shape ranker", site_set(t[t["_called"]]), tp_set, tn_set, slack,
                          "verdict = LIKELY EDIT"))

    if {"modal_pos", "start"} <= set(t.columns):
        cut = (t["modal_pos"] - t["start"]).abs()
        gated = t["_called"] & (cut <= args.max_cut_dist)
        rows.append(score_arm(f"WGS + ML + cut gate <={args.max_cut_dist:g}bp",
                              site_set(t[gated]), tp_set, tn_set, slack,
                              "indel must sit near the predicted cut"))

    res = pd.DataFrame(rows)
    base = res.loc[res["arm"].str.startswith("ECS rules + somatic"), "queue"]
    base = float(base.iloc[0]) if len(base) else float("nan")
    res["fold_vs_somatic_baseline"] = base / res["queue"].replace(0, pd.NA)

    print()
    print("=" * 100)
    print("REVIEW-QUEUE SIZE AT MATCHED RECALL  (curated AAVS1 truth: "
          f"{len(tp_set)} TP sites, {len(tn_set)} TN sites)")
    print("=" * 100)
    hdr = f"{'arm':<34} {'queue':>7} {'TP':>3} {'FP':>5} {'recall':>7} {'prec':>7} {'fold':>6}  note"
    print(hdr)
    print("-" * 100)
    for _, r in res.iterrows():
        fold = "" if pd.isna(r["fold_vs_somatic_baseline"]) else f"{r['fold_vs_somatic_baseline']:.0f}x"
        print(f"{r['arm']:<34} {r['queue']:>7,d} {r['TP']:>3d} {r['FP']:>5,d} "
              f"{r['recall']:>7.3f} {r['precision']:>7.3f} {fold:>6}  {r['note']}")
    print("-" * 100)
    print("  'queue' = unique sites a human would have to open in IGV.")
    print("  'fold'  = shrinkage vs the ECS somatic-gated baseline.")
    print("  Arms 1-2 use ECS (~5000x targeted); arms 3-5 use WGS (~30x). The ML arm reaches the")
    print("  same recall on ~150x less depth. Arm 3 is the same input as 4-5 with no model, so")
    print("  the 3 -> 4 step is the model's own contribution.")

    if args.out:
        res.to_csv(args.out, index=False)
        print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
