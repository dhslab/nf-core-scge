#!/usr/bin/env python3
"""
recall_vs_vaf.py — the honest limit, measured.

From the ECS⋈WGS training table, compute how often the WGS shape score recovers an
ECS-confirmed edit as a function of the ECS (error-corrected) VAF. This is THE
deliverable that keeps the WGS-only promise credible: it names the VAF above which
WGS-only detection is trustworthy, rather than implying WGS sees everything ECS sees.

WGS "detected" = the pipeline reported LIKELY EDIT (model score >= --hi, or the
high-evidence rescue). Sites with no spanning WGS reads are UNEVALUABLE rather than
missed, and are reported in their own column instead of being charged against recall —
the WGS depth floor is a property of the data, not of the scorer. `recall` is therefore
over evaluable sites; `recall_incl_unevaluable` keeps the pessimistic view.

THE DENOMINATOR (this is what makes the number mean anything)
------------------------------------------------------------
A recall figure is only as honest as the set of "real edits" it divides by. The raw ECS
label counts ANY nonzero indel fraction as an edit (`offtarget_ecs_edit_threshold`
defaults to 0.0), and at ECS depth that is overwhelmingly noise: in a real AAVS1 run,
10,780 of the 12,067 sites with any indel sat below 0.5% VAF with a MEDIAN of 3 indel
reads out of ~5,000. Dividing by those produced a headline recall of ~0.002 that said
nothing about whether the pipeline finds edits — it only measured how much ECS noise it
(correctly) ignores.

So the denominator here is a CREDIBLE ECS edit:
    label == 1                      somatic (ECS indel, absent from the matched normal)
    ecs_if          >= --min-ecs-vaf     above the ECS noise floor
    ecs_indel_reads >= --min-ecs-reads   with actual read support behind it
Read support is the load-bearing half: VAF alone cannot separate a genuine 0.5% edit at
5,000x from 3 stray reads, but read support can. Excluded sites are reported, and the
thresholds are written into the output CSV so the file is self-describing.
"""
import argparse
import numpy as np
import pandas as pd

BINS = [0.0, 0.005, 0.01, 0.02, 0.05, 0.10, 0.20, 0.50, 1.01]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--training", required=True, help="training.tsv from join_training_table.py")
    ap.add_argument("--hi", type=float, default=0.60, help="WGS score >= HI counts as detected")
    ap.add_argument("--detect-by-score", action="store_true",
                    help="score >= HI is the ONLY detection criterion; ignore the verdict "
                         "column (so high-evidence rescues do not count as detected)")
    ap.add_argument("--target-recall", type=float, default=0.80,
                    help="report the VAF floor where binned recall first reaches this")
    ap.add_argument("--min-ecs-vaf", type=float, default=0.005,
                    help="denominator: ECS VAF floor below which a call is assay noise, "
                         "not a real edit (0 = keep every ECS call)")
    ap.add_argument("--min-ecs-reads", type=int, default=5,
                    help="denominator: ECS indel reads required to believe the site is a "
                         "real edit rather than noise (0 = no read-support requirement)")
    ap.add_argument("--out-metrics", default="recall_vs_vaf.csv")
    ap.add_argument("--out-curve", default="recall_vs_vaf.png")
    args = ap.parse_args()

    df = pd.read_csv(args.training, sep="\t")
    somatic = df[df["label"] == 1].copy()

    # Restrict to CREDIBLE ECS edits (see module docstring). Anything failing these is
    # ECS noise, not a miss, and must not sit in the denominator.
    vaf = pd.to_numeric(somatic.get("ecs_if"), errors="coerce").fillna(0.0)
    keep = vaf >= args.min_ecs_vaf
    if "ecs_indel_reads" in somatic.columns:
        reads = pd.to_numeric(somatic["ecs_indel_reads"], errors="coerce").fillna(0)
        keep &= reads >= args.min_ecs_reads
    elif args.min_ecs_reads > 0:
        print(f"WARN: training table has no ecs_indel_reads column (produced before read "
              f"support was carried through) — the --min-ecs-reads {args.min_ecs_reads} "
              f"filter is INACTIVE and the denominator may still contain ECS noise.")
    pos = somatic[keep].copy()
    n_excluded = int(len(somatic) - len(pos))
    print(f"denominator: {len(pos)} credible ECS edits "
          f"(ecs_if>={args.min_ecs_vaf}, ecs_indel_reads>={args.min_ecs_reads}); "
          f"excluded {n_excluded} sub-threshold/no-support ECS calls as noise, "
          f"from {len(somatic)} somatic ECS-positive rows")

    if len(pos) == 0:
        print("no credible ECS edits in training table; nothing to measure")
        pd.DataFrame(columns=["vaf_bin", "n", "n_detected", "recall"]).to_csv(
            args.out_metrics, index=False)
        return

    pos["score"] = pd.to_numeric(pos["score"], errors="coerce")
    # "Detected" must mean what the pipeline actually REPORTS. The scorer can call a
    # LIKELY EDIT via the high-evidence rescue at a sub-HI model score, so scoring recall
    # on the raw model score alone under-counts exactly those recovered edits. Prefer the
    # verdict when it is present; fall back to the score for verdict-less tables.
    by_score = (pos["score"] >= args.hi).fillna(False)
    if "verdict" in pos.columns and not args.detect_by_score:
        pos["detected"] = pos["verdict"].astype(str).str.contains(
            "LIKELY EDIT", na=False).astype(int)
        n_extra = int((pos["detected"].astype(bool) & ~by_score).sum())
        if n_extra:
            print(f"detection = verdict contains 'LIKELY EDIT' "
                  f"({n_extra} rescued below score {args.hi})")
    else:
        pos["detected"] = by_score.astype(int)
    pos["vaf_bin"] = pd.cut(pd.to_numeric(pos["ecs_if"], errors="coerce"), bins=BINS,
                            right=False)

    # A site with no spanning WGS reads is UNEVALUABLE, not missed — the pipeline was
    # never given the chance to call it. Counting those as recall failures conflates the
    # WGS depth floor with scoring quality, so they are split out: `recall` is over sites
    # WGS could actually judge, and `recall_incl_unevaluable` keeps the pessimistic view.
    unevaluable = pos["verdict"].astype(str).str.contains(
        "INSUFFICIENT COVERAGE|NO CRAM", na=False) if "verdict" in pos.columns \
        else pd.Series(False, index=pos.index)
    pos["_uneval"] = unevaluable.astype(int)

    g = (pos.groupby("vaf_bin", observed=True)
            .agg(n=("detected", "size"), n_unevaluable=("_uneval", "sum"),
                 n_detected=("detected", "sum"))
            .reset_index())
    g["n_evaluable"] = g["n"] - g["n_unevaluable"]
    g["recall"] = (g["n_detected"] / g["n_evaluable"]).where(g["n_evaluable"] > 0)
    g["recall_incl_unevaluable"] = g["n_detected"] / g["n"]
    g = g[["vaf_bin", "n", "n_unevaluable", "n_evaluable", "n_detected",
           "recall", "recall_incl_unevaluable"]]
    # Stamp the denominator definition into the file itself. This CSV gets read months
    # later out of context; without these columns a reader cannot tell whether a low
    # recall means "missed real edits" or "divided by ECS noise".
    g["denom_min_ecs_vaf"] = args.min_ecs_vaf
    g["denom_min_ecs_reads"] = args.min_ecs_reads
    g["denom_excluded_as_noise"] = n_excluded
    g.to_csv(args.out_metrics, index=False)

    n_eval = int((~pos["_uneval"].astype(bool)).sum())
    overall = (pos.loc[~pos["_uneval"].astype(bool), "detected"].mean()
               if n_eval else float("nan"))
    floor = None
    for _, r in g.iterrows():
        if r["n"] >= 1 and r["recall"] >= args.target_recall:
            floor = r["vaf_bin"].left
            break
    print(f"overall WGS recall of credible ECS edits: {overall:.2f} over {n_eval} "
          f"EVALUABLE sites ({len(pos) - n_eval} of {len(pos)} had no WGS coverage)")
    print(f"VAF floor for >= {args.target_recall:.0%} binned recall: "
          f"{'>%.3f' % floor if floor is not None else 'not reached in these bins'}")
    print(g.to_string(index=False))

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        centers = [iv.left for iv in g["vaf_bin"]]
        fig, ax = plt.subplots(figsize=(7, 4.5))
        ax.plot(centers, g["recall"], "o-", color="#2b6cb0")
        ax.axhline(args.target_recall, ls="--", color="#a0aec0",
                   label=f"target recall {args.target_recall:.0%}")
        ax.set_xscale("symlog", linthresh=0.005)
        ax.set_xlabel("ECS error-corrected VAF (lower bin edge)")
        # NOT "score >= hi": detection is the reported verdict, which includes
        # high-evidence rescues below that score. And the ratio is over EVALUABLE sites.
        ax.set_ylabel("WGS recall (LIKELY EDIT / evaluable)")
        ax.set_ylim(-0.02, 1.02)
        # The figure travels further than the CSV (slides, papers), so it has to carry
        # its own denominator definition too.
        ax.set_title("WGS-only recovery of ECS-confirmed edits vs VAF")
        fig.text(0.5, 0.005,
                 f"denominator = {len(pos)} credible ECS edits "
                 f"(VAF ≥ {args.min_ecs_vaf:g}, ≥ {args.min_ecs_reads} ECS indel reads); "
                 f"{n_excluded} sub-threshold ECS calls excluded as assay noise",
                 ha="center", fontsize=7, color="#4a5568")
        for _, r in g.iterrows():
            # Annotate the EVALUABLE count — the denominator this point was actually
            # computed from. Labelling the credible total instead reads as "n/n detected"
            # and hides the depth floor (e.g. a bin of 4 credible edits with 2 uncovered
            # plots at recall 1.0 off 2 sites, not 4).
            lab = f"n={int(r['n_evaluable'])}"
            if r["n_unevaluable"]:
                lab += f" (+{int(r['n_unevaluable'])} uncov.)"
            ax.annotate(lab, (r["vaf_bin"].left, r["recall"]),
                        textcoords="offset points", xytext=(0, 6), fontsize=8, ha="center")
        ax.legend()
        fig.tight_layout()
        fig.savefig(args.out_curve, dpi=130)
        print(f"wrote {args.out_curve}")
    except Exception as e:  # plotting is a nicety; metrics CSV is the source of truth
        print(f"(curve not rendered: {e})")


if __name__ == "__main__":
    main()
