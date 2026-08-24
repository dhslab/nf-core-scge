#!/usr/bin/env python3
"""
train_shape_model.py — offline trainer for the off-target WGS shape model.

Reads a training.tsv (as emitted by BUILD_TRAINING_TABLE: WGS pileup features x
ECS-truth VAF + label) and fits the Stage-2 shape ranker consumed by score.py and
worklist_from_vcf.py. Emits a joblib bundle {model, features, ...} that is a drop-in
for --offtarget_shape_model.

The bundle records exactly which features the model was trained on; score.py and
worklist_from_vcf.py select model inputs by bundle["features"], so a training table
that carries only a subset of MODEL_FEATURES still yields a deployable (if weaker)
model. Patch score.py / join_training_table.py to carry the full feature set for the
strongest model.

Run inside the pipeline's off-target container so the pickle is written under the
same scikit-learn the pipeline runs — the version guard in features.check_sklearn_version
raises if the deployment runtime is OLDER than the pickle.
"""
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

# MODEL_FEATURES is defined once in features.py (the single source of truth). Import
# it if this script sits beside features.py (it does, in bin/); otherwise fall back to
# a literal copy kept in sync with features.MODEL_FEATURES.
try:
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from features import MODEL_FEATURES
except Exception:
    MODEL_FEATURES = ["indel_frac", "conc_ratio", "pos_conc", "pos_mad",
                      "modal_len", "modal_mapq", "softclip_frac"]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--training", required=True,
                    help="training.tsv from BUILD_TRAINING_TABLE (WGS features + label)")
    ap.add_argument("--out", default="wgs_shape_model.pkl", help="output model bundle (.pkl)")
    ap.add_argument("--metrics", default="train_metrics.json", help="output metrics JSON")
    ap.add_argument("--learning-rate", type=float, default=0.05)
    ap.add_argument("--max-iter", type=int, default=200)
    ap.add_argument("--max-depth", type=int, default=3)
    ap.add_argument("--max-leaf-nodes", type=int, default=31)
    ap.add_argument("--l2", type=float, default=1.0, help="L2 regularization")
    ap.add_argument("--holdout-frac", type=float, default=0.25,
                    help="stratified holdout fraction for reporting AUC/AP (0 = train on all, no holdout)")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--allow-nan-features", action="store_true",
                    help="keep rows with NaN features instead of dropping them. "
                         "HistGradientBoostingClassifier supports missing values natively; "
                         "needed whenever a feature is meaningfully absent (cut_dist is NaN "
                         "wherever no indel was observed, i.e. ~90%% of a hotspot panel).")
    args = ap.parse_args()

    import sklearn
    import joblib
    from sklearn.ensemble import HistGradientBoostingClassifier
    from sklearn.metrics import average_precision_score, roc_auc_score
    from sklearn.model_selection import train_test_split

    df = pd.read_csv(args.training, sep="\t")
    if "label" not in df.columns:
        sys.exit("ERROR: training table has no 'label' column")

    feats = [f for f in MODEL_FEATURES if f in df.columns]
    if not feats:
        sys.exit(f"ERROR: none of the model features {MODEL_FEATURES} are present in the "
                 f"training table (columns: {list(df.columns)})")
    missing = [f for f in MODEL_FEATURES if f not in df.columns]
    if missing:
        print(f"WARN: training table is missing {missing}; training a {len(feats)}-feature model "
              f"on {feats}. Patch score.py / join_training_table.py to carry the full feature set "
              f"for the strongest model.", file=sys.stderr)

    X = df[feats].apply(pd.to_numeric, errors="coerce")
    y = pd.to_numeric(df["label"], errors="coerce")
    # HistGradientBoostingClassifier handles NaN natively (it learns a missing-value branch
    # per split), so requiring every feature to be present is a choice, not a necessity --
    # and an expensive one now that the feature set includes cut_dist, which is legitimately
    # NaN at any site with no observed indel. On the CART panel that is 90% of rows
    # (9,278 of 99,308 have a cut_dist), so the default drops nine tenths of the training
    # data the moment cut_dist is added. --allow-nan-features keeps those rows and lets the
    # model treat "no indel to measure" as its own branch. Default is unchanged so existing
    # models stay reproducible.
    ok = y.notna()
    if not args.allow_nan_features:
        ok &= X.notna().all(axis=1)
    dropped = int((~ok).sum())
    X, y = X[ok].reset_index(drop=True), y[ok].astype(int).reset_index(drop=True)
    n_pos, n_neg = int((y == 1).sum()), int((y == 0).sum())
    if n_pos == 0 or n_neg == 0:
        sys.exit(f"ERROR: need both classes to train; got {n_pos} positive / {n_neg} negative rows")

    def make_clf():
        return HistGradientBoostingClassifier(
            learning_rate=args.learning_rate, max_iter=args.max_iter, max_depth=args.max_depth,
            max_leaf_nodes=args.max_leaf_nodes, l2_regularization=args.l2,
            random_state=args.seed, early_stopping=False)

    metrics = {"n_rows": int(len(y)), "n_pos": n_pos, "n_neg": n_neg, "n_dropped_incomplete": dropped,
               "features": feats, "missing_features": missing, "sklearn": sklearn.__version__,
               "hyperparams": {"learning_rate": args.learning_rate, "max_iter": args.max_iter,
                               "max_depth": args.max_depth, "max_leaf_nodes": args.max_leaf_nodes,
                               "l2_regularization": args.l2, "seed": args.seed}}

    # Optional stratified holdout to report generalization (AUC/AP) before the final fit.
    if 0 < args.holdout_frac < 1 and n_pos >= 2 and n_neg >= 2:
        X_tr, X_te, y_tr, y_te = train_test_split(
            X, y, test_size=args.holdout_frac, random_state=args.seed, stratify=y)
        p = make_clf().fit(X_tr, y_tr).predict_proba(X_te)[:, 1]
        metrics["holdout_n"] = int(len(y_te))
        metrics["holdout_auc"] = float(roc_auc_score(y_te, p)) if y_te.nunique() > 1 else None
        metrics["holdout_ap"] = float(average_precision_score(y_te, p))
    else:
        metrics["holdout_auc"] = None
        metrics["holdout_note"] = "holdout skipped (holdout_frac<=0 or too few of one class)"

    # Final model: fit on ALL rows for deployment.
    clf = make_clf().fit(X, y)

    bundle = {
        "model": clf,
        "features": feats,
        "role": "stage2_shape_ranker",
        "depth_augmented": False,
        "aug_depths": [],
        "stage1_filter": "is_target==0, min_mm<=2, indel_fraction>0.05, control_if<0.02, spanning>=8",
        "trained_on": (f"train_shape_model.py on {os.path.basename(args.training)} "
                       f"({n_pos} pos / {n_neg} neg loci, {len(feats)} features)"),
        "n_pos_loci": n_pos,
        "n_neg_loci": n_neg,
    }
    joblib.dump(bundle, args.out)
    with open(args.metrics, "w") as fh:
        json.dump(metrics, fh, indent=2)

    auc = metrics.get("holdout_auc")
    print(f"wrote {args.out}: HistGradientBoostingClassifier on {len(feats)} features "
          f"[{', '.join(feats)}], {n_pos} pos / {n_neg} neg"
          + (f", holdout AUC={auc:.3f}" if auc is not None else " (no holdout AUC)"))
    print(f"wrote {args.metrics}")


if __name__ == "__main__":
    main()
