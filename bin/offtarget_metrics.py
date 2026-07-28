#!/usr/bin/env python3
"""
offtarget_metrics.py — the metrics a recall-first diagnostic is actually judged on.

Reports PR-AUC and recall-weighted F-beta (F2, F5) for the WGS shape score, plus the
precision/recall/F1 at the operating point the pipeline actually reports (verdict =
LIKELY EDIT). ROC-AUC is carried for continuity but is deliberately NOT the headline:
at low positive prevalence it is dominated by the true-negative mass and flatters the
ranker. PR-AUC has no such property, which is why it leads here.

    F_beta = (1 + beta^2) * P * R / (beta^2 * P + R)

beta weights recall beta^2x more than precision (beta=2 -> 4x, beta=5 -> 25x). F_beta is
MONOTONE in beta: it rises with beta when R > P and falls when P > R, and always lies
between min(P, R) and max(P, R). So F1 is always an endpoint of {F1, F2, F5}, never the
middle value — if you see F1 in the middle, the beta wiring is inverted.

WHY THESE METRICS ARE NOT COMPUTED AGAINST THE MANUAL REVIEW
------------------------------------------------------------
This is the trap this script exists to avoid. The human-reviewed gold standard
(`cart_ecs_merged.csv.gz`) has 55 rows with manual_review == '1', 1 with '1?', and
80,440 NaN. **NaN means NOT REVIEWED, not reviewed-and-rejected.** There are no
confirmed negatives in it at all.

Precision therefore CANNOT be computed against that table. Treating NaN as negative
would count every genuine discovery the reviewers never got to as a false positive and
manufacture a confidently wrong number. That is why `validate_recall.py` reports recall
only, and why it stays the gold standard for recall. Do not "fix" it by filling NaN.

So there are two denominators, and this script keeps them strictly apart:

  1. RANKING + PRECISION  ->  training.tsv, which has a genuine two-class ECS label.
                              That is what this script computes. Precision here is
                              precision against the ECS label, NEVER against human review.
  2. RECALL vs HUMAN REVIEW -> `validate_recall.py`, reported as recall alone, unchanged.

Both denominators are stamped into the JSON and the text report so a reader months later
cannot mistake one for the other.

THE POSITIVE SET (same credibility gates as the recall curve)
-------------------------------------------------------------
A raw `label == 1` counts ANY nonzero ECS indel fraction as an edit, which at ECS depth
is overwhelmingly noise. Positives are therefore gated exactly as `recall_vs_vaf.py`
gates its denominator:

    label == 1  AND  ecs_if >= --min-ecs-vaf  AND  ecs_indel_reads >= --min-ecs-reads

`label == 1` rows failing those gates are AMBIGUOUS, not negative: ECS saw something,
below the credibility floor. They are EXCLUDED and counted, for the same reason the
unreviewed manual-review rows are excluded — asserting they are negatives is an
assumption the data does not support. `--ambiguous-as-negative` reports the alternative,
and the JSON always carries it as a labelled sensitivity block so the choice is visible.

THE NEGATIVE SET
----------------
`label == 0` is two different things:
  * ecs_is_edit == 0  — a genuine ECS-negative at a scored hotspot (a true negative);
  * ecs_is_edit == 1  — an ECS edit DEMOTED because the indel is in the matched normal
                        (germline/artifact).
The germline-demoted rows carry a real indel, so the shape model correctly scores their
pileup shape high — germline rejection is a separate downstream gate (their verdict is
`GERMLINE/ARTIFACT (in normal)`, and none are called LIKELY EDIT). Ranking the shape
score against them charges it with a job it does not do and is not asked to do, so the
default negative set is `ecs_negative`. `--negatives all_label0` gives the stricter view
and is reported as a sensitivity block either way.
"""
import argparse
import json
import sys

import numpy as np
import pandas as pd

# VAF floors for the stratified view. A single PR-AUC over all credible positives is
# dominated by sub-1% VAF sites that WGS at ~30x genuinely cannot see, so it understates
# the ranker exactly where the ranker is usable. Stratifying by ECS VAF separates
# "the model cannot rank" from "the data has no signal at this depth".
VAF_FLOORS = [0.005, 0.01, 0.02, 0.05, 0.10]

DETECT_PATTERN = "LIKELY EDIT"
UNEVALUABLE_PATTERN = "INSUFFICIENT COVERAGE|NO CRAM"

MANUAL_REVIEW_NOTE = (
    "Precision is NOT computable against the human manual review: that table has 55 "
    "confirmed positives ('1'), 1 '1?', and 80,440 NaN, where NaN means UNREVIEWED, not "
    "rejected. It contains no confirmed negatives, so any precision computed from it "
    "would score every discovery the reviewers never reached as a false positive. All "
    "precision/PR-AUC/F-beta figures in this file use the ECS label in training.tsv as "
    "their denominator. Recall against human review is reported separately, and as "
    "recall only, by bin/validate_recall.py."
)


def fbeta(precision, recall, beta):
    """F_beta = (1 + b^2) * P * R / (b^2 * P + R); 0 when the denominator vanishes."""
    denom = (beta * beta * precision) + recall
    if denom <= 0:
        return 0.0
    return float((1.0 + beta * beta) * precision * recall / denom)


def build_masks(df, args):
    """Split the training table into credible positives / negatives / ambiguous."""
    label = pd.to_numeric(df.get("label"), errors="coerce")
    vaf = pd.to_numeric(df.get("ecs_if"), errors="coerce").fillna(0.0)

    gate = vaf >= args.min_ecs_vaf
    warn = None
    if "ecs_indel_reads" in df.columns:
        reads = pd.to_numeric(df["ecs_indel_reads"], errors="coerce").fillna(0)
        gate &= reads >= args.min_ecs_reads
    elif args.min_ecs_reads > 0:
        # Same guard as recall_vs_vaf.py: a table built before read support was carried
        # through silently disables half the credibility gate. Warn, never pass quietly.
        warn = (f"training table has no ecs_indel_reads column (produced before read "
                f"support was carried through) — the --min-ecs-reads {args.min_ecs_reads} "
                f"filter is INACTIVE and the positive set may still contain ECS noise")
        print(f"WARN: {warn}", file=sys.stderr)

    pos = (label == 1) & gate
    ambiguous = (label == 1) & ~gate

    if args.negatives == "ecs_negative" and "ecs_is_edit" in df.columns:
        ecs_edit = pd.to_numeric(df["ecs_is_edit"], errors="coerce").fillna(0)
        neg = (label == 0) & (ecs_edit == 0)
        neg_germline = (label == 0) & (ecs_edit != 0)
    else:
        neg = (label == 0)
        neg_germline = pd.Series(False, index=df.index)

    if args.ambiguous_as_negative:
        neg = neg | ambiguous
        ambiguous = pd.Series(False, index=df.index)

    return pos, neg, ambiguous, neg_germline, warn


def score_and_calls(df, hi):
    """Continuous score + the reported operating point + the unevaluable mask.

    The operating point must be what the pipeline actually REPORTS. The scorer can call
    a LIKELY EDIT via the high-evidence rescue at a sub-`hi` model score, so scoring the
    raw threshold alone would under-count exactly those recovered edits. Prefer the
    verdict; fall back to `score >= hi` only for verdict-less tables. Same precedence as
    recall_vs_vaf.py, so the two files agree on what "detected" means.
    """
    score = pd.to_numeric(df.get("score"), errors="coerce")
    if "verdict" in df.columns:
        verdict = df["verdict"].astype(str)
        called = verdict.str.contains(DETECT_PATTERN, na=False)
        uneval = verdict.str.contains(UNEVALUABLE_PATTERN, na=False)
        criterion = f"verdict contains '{DETECT_PATTERN}'"
    else:
        called = (score >= hi).fillna(False)
        uneval = score.isna()
        criterion = f"score >= {hi} (table has no verdict column)"
    return score, called, uneval, criterion


def ranking_metrics(y, s):
    """PR-AUC / ROC-AUC on a finite-score, two-class subset. None when undefined."""
    from sklearn.metrics import average_precision_score, roc_auc_score
    out = {"pr_auc": None, "roc_auc": None}
    if len(y) == 0 or len(np.unique(y)) < 2:
        return out
    out["pr_auc"] = float(average_precision_score(y, s))
    out["roc_auc"] = float(roc_auc_score(y, s))
    return out


def evaluate(df, pos, neg, score, called, uneval, betas, criterion=None):
    """Full metric block for one positive/negative definition."""
    m = pos | neg
    y_all = pos[m].astype(int).to_numpy()
    s_all = score[m].to_numpy()
    called_all = called[m].to_numpy()
    uneval_all = uneval[m].to_numpy()

    # Rows with no score (INSUFFICIENT COVERAGE / NO CRAM) are excluded from the ranking
    # metrics: average_precision_score needs a continuous score, and filling them with 0
    # would fabricate confident negatives. They are counted, not hidden.
    finite = np.isfinite(s_all)
    y, s, pred = y_all[finite], s_all[finite], called_all[finite].astype(int)

    n_pos, n_neg = int((y == 1).sum()), int((y == 0).sum())
    block = {
        "n_pos": n_pos,
        "n_neg": n_neg,
        "n_total": n_pos + n_neg,
        "prevalence": float(y.mean()) if len(y) else None,
        "n_excluded_no_score": int((~finite).sum()),
        "n_excluded_no_score_pos": int((y_all[~finite] == 1).sum()),
        "n_excluded_no_score_neg": int((y_all[~finite] == 0).sum()),
        "n_unevaluable": int(uneval_all.sum()),
    }
    block.update(ranking_metrics(y, s))

    tp = int(((pred == 1) & (y == 1)).sum())
    fp = int(((pred == 1) & (y == 0)).sum())
    fn = int(((pred == 0) & (y == 1)).sum())
    tn = int(((pred == 0) & (y == 0)).sum())
    precision = tp / (tp + fp) if (tp + fp) else 0.0
    recall = tp / (tp + fn) if (tp + fn) else 0.0

    op = {
        "criterion": criterion or f"verdict contains '{DETECT_PATTERN}'",
        "tp": tp, "fp": fp, "fn": fn, "tn": tn,
        "precision": float(precision),
        "recall": float(recall),
        "f1": fbeta(precision, recall, 1.0),
        # Pessimistic view: charge the no-score positives against recall as well, the
        # same split recall_vs_vaf.csv makes with recall vs recall_incl_unevaluable.
        "recall_incl_unevaluable": (
            float(tp / int((y_all == 1).sum())) if int((y_all == 1).sum()) else None),
    }
    for b in betas:
        op[f"f{_beta_key(b)}"] = fbeta(precision, recall, b)
    block["operating_point"] = op

    # PR-AUC below prevalence is the "useless ranker or a bug" tripwire from the design
    # notes. Here it is usually neither: it is the WGS depth floor showing up, because
    # the credible-positive set is dominated by sub-1% VAF sites. Say so rather than
    # emit a bare number that reads as a broken model.
    if block["pr_auc"] is not None and block["prevalence"] is not None:
        block["pr_auc_lift_over_prevalence"] = (
            float(block["pr_auc"] / block["prevalence"]) if block["prevalence"] else None)
    return block


def _beta_key(b):
    """2.0 -> '2', 2.5 -> '2.5' — keeps the common case reading as f2/f5."""
    return str(int(b)) if float(b).is_integer() else str(b)


def vaf_stratified(df, pos, neg, score, called, uneval, betas, args):
    """PR-AUC vs prevalence as the ECS VAF floor on the POSITIVE set rises."""
    vaf = pd.to_numeric(df.get("ecs_if"), errors="coerce").fillna(0.0)
    rows = []
    for floor in VAF_FLOORS:
        if floor < args.min_ecs_vaf:
            continue
        p = pos & (vaf >= floor)
        if int(p.sum()) == 0:
            continue
        b = evaluate(df, p, neg, score, called, uneval, betas)
        rows.append({
            "ecs_vaf_floor": floor,
            "n_pos": b["n_pos"], "n_neg": b["n_neg"],
            "prevalence": b["prevalence"], "pr_auc": b["pr_auc"],
            "pr_auc_lift_over_prevalence": b.get("pr_auc_lift_over_prevalence"),
            "roc_auc": b["roc_auc"],
            "recall": b["operating_point"]["recall"],
            "precision": b["operating_point"]["precision"],
        })
    return rows


def fmt(v, spec=".4f"):
    return "n/a" if v is None or (isinstance(v, float) and not np.isfinite(v)) \
        else format(v, spec)


def render_text(out, betas):
    """One screen. Every number carries the denominator it was computed against."""
    m = out["metrics"]
    op = m["operating_point"] if m else None
    L = []
    A = L.append
    A("=" * 78)
    A("OFF-TARGET CLINICAL METRICS — PR-AUC and recall-weighted F-beta")
    A("=" * 78)
    A(f"training table : {out['inputs']['training']}")
    A(f"rows           : {out['inputs']['n_rows']}")
    A("")
    A("-- DENOMINATOR 1 of 2: the ECS label (this is what every number below uses) --")
    A(f"  positives : label==1 AND ecs_if >= {out['inputs']['min_ecs_vaf']} "
      f"AND ecs_indel_reads >= {out['inputs']['min_ecs_reads']}")
    A(f"  negatives : {out['inputs']['negatives']}"
      f"{'  (label==0 AND ecs_is_edit==0)' if out['inputs']['negatives'] == 'ecs_negative' else '  (all label==0)'}")
    counts = out.get("counts") or {}
    A(f"  ambiguous : {counts.get('n_ambiguous', 0)} label==1 rows below the credibility "
      f"floor — EXCLUDED, not counted as negatives")
    if counts.get("n_germline_demoted_excluded"):
        A(f"  excluded  : {counts['n_germline_demoted_excluded']} germline-demoted "
          f"label==0 rows (indel in the matched normal; rejected by a separate gate)")
    A("")
    A("-- DENOMINATOR 2 of 2: the human manual review --")
    for line in _wrap(MANUAL_REVIEW_NOTE, 74):
        A(f"  {line}")
    A("")
    if not m:
        A(f"NOTE: {out['note']}")
        A("=" * 78)
        return "\n".join(L) + "\n"

    A("-- RANKING (score as a continuous ranker; ECS-label denominator) --")
    A(f"  n_pos / n_neg      : {m['n_pos']} / {m['n_neg']}   prevalence {fmt(m['prevalence'])}")
    A(f"  PR-AUC             : {fmt(m['pr_auc'])}"
      + (f"   ({fmt(m.get('pr_auc_lift_over_prevalence'), '.2f')}x prevalence)"
         if m.get("pr_auc_lift_over_prevalence") is not None else ""))
    A(f"  ROC-AUC            : {fmt(m['roc_auc'])}   (kept for continuity; not the headline "
      f"at low prevalence)")
    A(f"  excluded, no score : {m['n_excluded_no_score']} "
      f"({m['n_excluded_no_score_pos']} pos / {m['n_excluded_no_score_neg']} neg) — "
      f"INSUFFICIENT COVERAGE / NO CRAM, never filled with 0")
    A("")
    A(f"-- OPERATING POINT ({op['criterion']}; ECS-label denominator) --")
    A(f"  TP {op['tp']}   FP {op['fp']}   FN {op['fn']}   TN {op['tn']}")
    A(f"  precision (vs ECS label, NOT vs human review) : {fmt(op['precision'])}")
    A(f"  recall                                        : {fmt(op['recall'])}")
    A(f"  recall incl. unevaluable                      : {fmt(op['recall_incl_unevaluable'])}")
    A(f"  {'F1':<44s}: {fmt(op['f1'])}")
    for b in betas:
        k = f"f{_beta_key(b)}"
        A(f"  {f'F{_beta_key(b)} (recall weighted {b * b:g}x)':<44s}: {fmt(op[k])}")
    A("")
    if out["by_vaf_floor"]:
        A("-- BY ECS VAF FLOOR ON THE POSITIVE SET --")
        A("   (a single PR-AUC over all credible positives is dominated by sub-1% VAF")
        A("    sites that WGS at ~30x cannot see; this separates ranker from depth floor)")
        A(f"   {'vaf>=':>7} {'n_pos':>6} {'n_neg':>6} {'prev':>7} {'PR-AUC':>7} {'lift':>6} {'ROC':>6}")
        for r in out["by_vaf_floor"]:
            A(f"   {r['ecs_vaf_floor']:>7g} {r['n_pos']:>6} {r['n_neg']:>6} "
              f"{fmt(r['prevalence'], '.3f'):>7} {fmt(r['pr_auc'], '.3f'):>7} "
              f"{fmt(r.get('pr_auc_lift_over_prevalence'), '.2f'):>6} "
              f"{fmt(r['roc_auc'], '.3f'):>6}")
        A("")
    if out["sensitivity"]:
        A("-- SENSITIVITY: the same score under other defensible denominators --")
        for s in out["sensitivity"]:
            A(f"   {s['definition']:<34} n_pos {s['n_pos']:>5}  n_neg {s['n_neg']:>5}  "
              f"prev {fmt(s['prevalence'], '.3f')}  PR-AUC {fmt(s['pr_auc'], '.3f')}")
        A("")
    if out["notes"]:
        A("-- NOTES --")
        for n in out["notes"]:
            for line in _wrap(n, 74):
                A(f"  {line}")
        A("")
    A("=" * 78)
    return "\n".join(L) + "\n"


def _wrap(text, width):
    words, line, lines = text.split(), "", []
    for w in words:
        if line and len(line) + 1 + len(w) > width:
            lines.append(line)
            line = w
        else:
            line = f"{line} {w}".strip()
    if line:
        lines.append(line)
    return lines


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--training", required=True,
                    help="training.tsv from join_training_table.py (needs label + score)")
    ap.add_argument("--hi", type=float, default=0.60,
                    help="score threshold; used only when the table has no verdict column")
    ap.add_argument("--min-ecs-vaf", type=float, default=0.005,
                    help="positive set: ECS VAF floor below which a call is assay noise")
    ap.add_argument("--min-ecs-reads", type=int, default=5,
                    help="positive set: ECS indel reads required to believe a site")
    ap.add_argument("--beta", type=float, action="append", default=None, metavar="B",
                    help="F-beta weight; repeatable (default: 2 and 5)")
    ap.add_argument("--betas", help="comma-separated F-beta weights, e.g. '2,5'")
    ap.add_argument("--negatives", choices=["ecs_negative", "all_label0"],
                    default="ecs_negative",
                    help="ecs_negative (default) = label==0 AND ecs_is_edit==0; "
                         "all_label0 also includes germline-demoted ECS edits")
    ap.add_argument("--ambiguous-as-negative", action="store_true",
                    help="treat sub-credibility label==1 rows as negatives (NOT the "
                         "default: ECS saw something there, so calling them negative is "
                         "an assumption the data does not support)")
    ap.add_argument("--out-json", default="offtarget_metrics.json")
    ap.add_argument("--out-txt", default="offtarget_metrics.txt")
    ap.add_argument("--out-curve", default=None, help="optional PR-curve PNG")
    args = ap.parse_args()

    betas = []
    if args.betas:
        betas += [float(b) for b in args.betas.replace(" ", "").split(",") if b]
    if args.beta:
        betas += list(args.beta)
    betas = sorted(set(betas)) or [2.0, 5.0]

    df = pd.read_csv(args.training, sep="\t")
    notes = []

    out = {
        "inputs": {
            "training": args.training,
            "n_rows": int(len(df)),
            "min_ecs_vaf": args.min_ecs_vaf,
            "min_ecs_reads": args.min_ecs_reads,
            "negatives": args.negatives,
            "ambiguous_as_negative": bool(args.ambiguous_as_negative),
            "betas": betas,
            "hi": args.hi,
        },
        "denominators": {
            "metrics_denominator": (
                "ECS label in training.tsv: positives are label==1 gated to credible ECS "
                "edits (ecs_if >= min_ecs_vaf AND ecs_indel_reads >= min_ecs_reads); "
                f"negatives are '{args.negatives}'."),
            "manual_review_denominator": MANUAL_REVIEW_NOTE,
        },
        "counts": {}, "metrics": None, "by_vaf_floor": [], "sensitivity": [],
        "notes": notes, "note": None,
    }

    if "label" not in df.columns or "score" not in df.columns:
        out["note"] = ("training table lacks a 'label' and/or 'score' column — no metrics "
                       "computable")
        _emit(out, betas, args)
        return

    pos, neg, ambiguous, neg_germline, warn = build_masks(df, args)
    if warn:
        notes.append(warn)
    score, called, uneval, criterion = score_and_calls(df, args.hi)

    out["counts"] = {
        "n_positives_credible": int(pos.sum()),
        "n_negatives": int(neg.sum()),
        "n_ambiguous": int(ambiguous.sum()),
        "n_germline_demoted_excluded": int(neg_germline.sum()),
        "n_somatic_label1_total": int((pd.to_numeric(df["label"], errors="coerce") == 1).sum()),
    }

    # Degenerate cases are legitimate results, not errors: a run with no credible edits
    # still has to publish a file. Emit nulls plus an explanatory note and exit 0.
    if int(pos.sum()) == 0 or int(neg.sum()) == 0:
        out["note"] = (
            f"single-class evaluation set ({int(pos.sum())} credible positives, "
            f"{int(neg.sum())} negatives) — PR-AUC and F-beta are undefined. This is a "
            f"legitimate outcome for a run with no credible ECS edits, not a failure.")
        _emit(out, betas, args)
        return

    out["metrics"] = evaluate(df, pos, neg, score, called, uneval, betas, criterion)

    if out["metrics"]["n_pos"] == 0 or out["metrics"]["n_neg"] == 0:
        out["note"] = ("every row of one class had a NaN score (INSUFFICIENT COVERAGE / "
                       "NO CRAM) — ranking metrics undefined on the evaluable subset")

    out["by_vaf_floor"] = vaf_stratified(df, pos, neg, score, called, uneval, betas, args)

    # Sensitivity: show what the other defensible denominators would give, so the one
    # judgement call in here is visible rather than buried in a default.
    for definition, p2, n2 in _alternatives(df, pos, neg, ambiguous, neg_germline, args):
        b = evaluate(df, p2, n2, score, called, uneval, betas)
        out["sensitivity"].append({
            "definition": definition, "n_pos": b["n_pos"], "n_neg": b["n_neg"],
            "prevalence": b["prevalence"], "pr_auc": b["pr_auc"], "roc_auc": b["roc_auc"],
        })

    m = out["metrics"]
    if m["pr_auc"] is not None and m["prevalence"] and m["pr_auc"] <= m["prevalence"]:
        notes.append(
            f"PR-AUC ({m['pr_auc']:.3f}) does not exceed prevalence "
            f"({m['prevalence']:.3f}): over this positive set the score adds no ranking "
            f"power. Check the by-VAF-floor table before concluding the model is broken "
            f"— the credible-positive set is usually dominated by sub-1% VAF sites that "
            f"WGS at ~30x cannot see, and the lift typically rises sharply with the VAF "
            f"floor. A flat lift at every floor does indicate a real problem.")

    _emit(out, betas, args)


def _alternatives(df, pos, neg, ambiguous, neg_germline, args):
    """(label, positives, negatives) triples for the sensitivity block."""
    alts = []
    if int(neg_germline.sum()):
        alts.append(("negatives = all label==0", pos, neg | neg_germline))
    if int(ambiguous.sum()) and not args.ambiguous_as_negative:
        alts.append(("+ sub-credibility as negative", pos, neg | ambiguous))
    label = pd.to_numeric(df.get("label"), errors="coerce")
    ungated = (label == 1)
    if int(ungated.sum()) != int(pos.sum()):
        alts.append(("ungated label==1 (ECS noise incl.)", ungated, neg))
    return alts


def _emit(out, betas, args):
    with open(args.out_json, "w") as fh:
        json.dump(out, fh, indent=2)
    text = render_text(out, betas)
    with open(args.out_txt, "w") as fh:
        fh.write(text)
    print(text)
    print(f"wrote {args.out_json} and {args.out_txt}")
    if args.out_curve:
        _plot(out, args)


def _plot(out, args):
    """PR curve with the operating point marked. A nicety — the JSON is the truth."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from sklearn.metrics import precision_recall_curve

        df = pd.read_csv(args.training, sep="\t")
        pos, neg, _amb, _gd, _w = build_masks(df, args)
        score, called, _u, _c = score_and_calls(df, args.hi)
        m = pos | neg
        y = pos[m].astype(int).to_numpy()
        s = score[m].to_numpy()
        ok = np.isfinite(s)
        y, s, pred = y[ok], s[ok], called[m].to_numpy()[ok]
        if len(np.unique(y)) < 2:
            print("(PR curve not rendered: single-class evaluation set)")
            return
        p, r, _ = precision_recall_curve(y, s)
        op = out["metrics"]["operating_point"]
        fig, ax = plt.subplots(figsize=(7, 4.5))
        ax.plot(r, p, "-", color="#2b6cb0", label=f"PR-AUC = {out['metrics']['pr_auc']:.3f}")
        ax.axhline(out["metrics"]["prevalence"], ls="--", color="#a0aec0",
                   label=f"prevalence = {out['metrics']['prevalence']:.3f}")
        ax.plot([op["recall"]], [op["precision"]], "o", color="#c53030", ms=9,
                label=f"LIKELY EDIT (P={op['precision']:.2f}, R={op['recall']:.2f})")
        ax.set_xlabel("recall")
        ax.set_ylabel("precision (vs ECS label — NOT vs human review)")
        ax.set_xlim(-0.02, 1.02)
        ax.set_ylim(-0.02, 1.02)
        ax.set_title("WGS shape score — precision/recall vs the ECS label")
        fig.text(0.5, 0.005,
                 f"denominator = {out['metrics']['n_pos']} credible ECS edits vs "
                 f"{out['metrics']['n_neg']} ECS-negatives; precision here is against the "
                 f"ECS label, never against human review",
                 ha="center", fontsize=7, color="#4a5568")
        ax.legend(loc="upper right", fontsize=8)
        fig.tight_layout()
        fig.savefig(args.out_curve, dpi=130)
        print(f"wrote {args.out_curve}")
    except Exception as e:
        print(f"(PR curve not rendered: {e})")


if __name__ == "__main__":
    main()
