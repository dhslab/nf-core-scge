#!/usr/bin/env python3
"""
join_training_table.py — the ECS⋈WGS join that makes the WGS-only model possible.

Left  : WGS per-hotspot features + shape score (score.py output; sample = WGS sample).
Right : ECS per-(guide,site) truth VAF + label (from hotspot_to_table.py).
Key   : (guide, chrom, start). WGS sample -> guide via the samplesheet.

Output: training.tsv, one row per (guide, WGS sample, hotspot) = WGS pileup features
(indel_frac, conc_ratio, score, spanning, ctrl_if, modal_len, verdict) + ecs_if (truth
VAF) + label. This feeds the OFFLINE trainer; the deployed model stays a fixed asset.
"""
import sys, argparse
import pandas as pd

WGS_FEATURE_COLS = ["indel_frac", "conc_ratio", "spanning", "ctrl_if", "modal_len",
                    "modal_pos", "min_mm", "score", "verdict"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--wgs-scores", required=True, help="score.py output CSV (WGS at hotspots)")
    ap.add_argument("--truth", required=True, help="ecs_hotspot_truth.csv from hotspot_to_table.py")
    ap.add_argument("--samplesheet", required=True)
    ap.add_argument("--out", default="training.tsv")
    args = ap.parse_args()

    ss = pd.read_csv(args.samplesheet)
    ss.columns = [c.strip().lower() for c in ss.columns]
    sample_guide = dict(zip(ss["sample"], ss["guide"]))

    wgs = pd.read_csv(args.wgs_scores)
    wgs["guide"] = wgs["sample"].map(sample_guide)
    keep = ["sample", "guide", "chrom", "start"] + [c for c in WGS_FEATURE_COLS if c in wgs.columns]
    wgs = wgs[keep].copy()

    truth = pd.read_csv(args.truth)  # guide, chrom, start, end, is_target, target_info, ecs_if, ecs_is_edit

    # normalize key dtypes
    for d in (wgs, truth):
        d["chrom"] = d["chrom"].astype(str)
        d["start"] = pd.to_numeric(d["start"], errors="coerce").astype("Int64")

    merged = wgs.merge(truth[["guide", "chrom", "start", "ecs_if", "ecs_is_edit", "is_target"]],
                       on=["guide", "chrom", "start"], how="inner")
    merged = merged.rename(columns={"ecs_is_edit": "label"})
    merged.to_csv(args.out, sep="\t", index=False)

    n_pos = int((merged["label"] == 1).sum())
    print(f"wrote {args.out}: {len(merged)} rows "
          f"({n_pos} ECS-positive, {len(merged) - n_pos} ECS-negative) "
          f"across {merged['guide'].nunique()} guides")
    if len(merged) == 0:
        msg = ("empty join — WGS sample->guide and hotspot coords do not line up with the "
               "ECS truth (a coordinate off-by-one between the arms is the usual culprit)")
        # If BOTH arms produced rows but nothing joined, this is a keying bug, not a
        # legitimately empty run: fail loudly so a silent empty training.tsv can't pass.
        # Only warn when one side is genuinely empty (e.g. an ECS- or WGS-less run).
        if len(wgs) > 0 and len(truth) > 0:
            sys.exit(f"ERROR: {msg} [{len(wgs)} WGS rows x {len(truth)} ECS truth rows -> 0 joined]")
        print(f"WARN: {msg} [{len(wgs)} WGS rows, {len(truth)} ECS truth rows]", file=sys.stderr)


if __name__ == "__main__":
    main()
