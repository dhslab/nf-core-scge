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

WGS_FEATURE_COLS = ["indel_frac", "conc_ratio", "pos_conc", "pos_mad", "modal_len",
                    "modal_mapq", "softclip_frac", "spanning", "ctrl_if",
                    "modal_pos", "min_mm", "score", "verdict", "call_basis"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--wgs-scores", required=True, help="score.py output CSV (WGS at hotspots)")
    ap.add_argument("--truth", required=True, help="ecs_hotspot_truth.csv from hotspot_to_table.py")
    ap.add_argument("--samplesheet", required=True)
    ap.add_argument("--germline-max-ctrl-if", type=float, default=0.05,
                    help="a hotspot whose matched-normal (WGS) indel fraction exceeds this is "
                         "germline/artifact, not a somatic edit -> label 0 even if ECS shows an indel")
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

    truth_cols = ["guide", "chrom", "start", "ecs_if", "ecs_is_edit", "is_target"]
    # read support is optional (older ecs_hotspot_truth.csv predates it) but is what
    # lets recall_vs_vaf.py keep ECS noise out of its denominator
    truth_cols += [c for c in ("ecs_indel_reads", "ecs_total_reads") if c in truth.columns]
    merged = wgs.merge(truth[truth_cols], on=["guide", "chrom", "start"], how="inner")

    # Training label = a SOMATIC edit: ECS saw an indel AND it is absent from the matched WGS
    # normal. ECS alone can't exclude germline here — its control is a different individual, so
    # control_indel_fraction is ~0 everywhere — so the matched-normal WGS signal (ctrl_if) is the
    # only thing that separates a real edit from a germline/recurrent variant. The raw ECS call is
    # kept as `ecs_is_edit`; `label` (what the trainer + recall metric consume) is the somatic one.
    is_edit = merged["ecs_is_edit"] == 1
    if "ctrl_if" in merged.columns:
        ctrl = pd.to_numeric(merged["ctrl_if"], errors="coerce").fillna(0.0)
        germline = is_edit & (ctrl > args.germline_max_ctrl_if)
        merged["label"] = (is_edit & (ctrl <= args.germline_max_ctrl_if)).astype(int)
    else:
        print("WARN: no ctrl_if (matched-normal) column — cannot exclude germline; "
              "label = raw ECS edit", file=sys.stderr)
        germline = pd.Series(False, index=merged.index)
        merged["label"] = is_edit.astype(int)
    merged.to_csv(args.out, sep="\t", index=False)

    n_pos = int((merged["label"] == 1).sum())
    n_germ = int(germline.sum())
    print(f"wrote {args.out}: {len(merged)} rows "
          f"({n_pos} somatic-edit positives, {n_germ} ECS edits demoted as germline/in-normal, "
          f"{len(merged) - n_pos - n_germ} ECS-negative) across "
          f"{merged['guide'].nunique() if len(merged) else 0} guides")
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
