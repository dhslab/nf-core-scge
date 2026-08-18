#!/usr/bin/env python3
"""cut_distance_model.py — fit the distance-from-cut distribution empirically, as a spatial prior.

The question
------------
`review_filter.py` treats distance from the predicted cut as a boolean: keep if `cut_dist_min <=
10`, drop otherwise. That threshold was never derived from the data — it is inherited from the
caller's own `-d/--max-mutation-distance` default. This script asks what the distribution actually
looks like and what a properly-specified spatial model would say instead.

The model
---------
Two components, one for each physical process:

  P(d | background)  an artifact has no reason to prefer any offset relative to a predicted cut
                     site, so the null is Uniform over the admissible window: P(d) = 1/W.
  P(d | signal)      Cas9 cuts blunt, 3 bp 5' of the PAM. A real edit therefore sits AT the cut,
                     broadened only by repair microhomology and by where the aligner chose to
                     place the indel. The prediction is a sharp spike at 0-1 bp.

Given both, the evidence carried by an observed distance is the likelihood ratio

    LR(d) = P(d | signal) / P(d | background)

which is a real likelihood ratio, unlike AQ -- see docs/NOISE_MODEL.md, which is careful about
that distinction. LR turns a boolean gate into a continuous quantity on the same evidence scale as
AQ, so the two can be added in log space instead of acting as independent hurdles.

Estimation
----------
The signal pmf is the empirical distribution with add-one smoothing (so no bin has zero
probability and LR stays finite). A discretised half-normal is also fitted and reported, because a
one-parameter form is easier to defend and to transfer to another cohort than 11 empirical bins.
The background is Uniform by construction; the script tests that assumption rather than asserting
it, with a chi-square over the observed off-target bins.

TRUNCATION — read this before quoting anything about the tail
--------------------------------------------------------------
The caller's `-d` default is 10 and `get_indels.nf` passes no override, so `min_cut_distance` is
only ever observed on events the caller already accepted at <= 10 bp. The window is 0..10 and the
distribution beyond it is not estimable from these tables: what is measured is the shape INSIDE
the window we are allowed to see, not the full spatial distribution.

`indel_info` field 7 is a different and always-larger quantity -- the distance from the event's
ANCHOR rather than the minimum over the event's span -- and it is not capped, reaching 28 bp. It
is reported here for tail shape only, clearly separated, and must not be pooled with
`min_cut_distance`.

usage:
  cut_distance_model.py --tables 'results_cart_bnd/*/*.offtarget_analysis.tsv' \
      --queue-all results_cart_bnd/review/review_queue_all.tsv \
      --truth-wgs '.../cart_wgs_merged.xlsx' --figdir docs/images
"""
import argparse
import glob
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import chisquare, norm

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import noise_model as nm                                              # noqa: E402

WINDOW = 11          # admissible bins, 0..10 inclusive -- set by the caller's -d default


def hist(d, w=WINDOW):
    d = pd.to_numeric(pd.Series(d), errors="coerce").dropna().astype(int)
    d = d[(d >= 0) & (d < w)]
    return np.array([int((d == i).sum()) for i in range(w)])


def smoothed_pmf(counts):
    """Add-one (Laplace) smoothing. Keeps LR finite in bins no signal event happened to land in."""
    c = counts.astype(float) + 1.0
    return c / c.sum()


def halfnormal_pmf(sigma, w=WINDOW):
    """Discretised half-normal on 0..w-1: mass of [d-0.5, d+0.5) folded about zero."""
    edges = np.arange(w + 1) - 0.5
    edges[0] = 0.0
    cdf = 2 * norm.cdf(edges / sigma) - 1.0
    p = np.diff(cdf)
    return p / p.sum()


def fit_halfnormal(d):
    """MLE for a half-normal scale is sqrt(mean(d^2)); the continuity correction is left out
    deliberately, since it moves sigma by less than the sampling error at these counts."""
    d = np.asarray(d, dtype=float)
    return float(np.sqrt(np.mean(d ** 2))) if len(d) else float("nan")


def sweep(pos, neg, label_pos="confirmed", label_neg="rejected"):
    print(f"    threshold   {label_pos} kept        {label_neg} kept")
    for t in range(0, WINDOW):
        p = int((pos <= t).sum())
        n = int((neg <= t).sum())
        print(f"      d<={t:<2d}      {p:3d}/{len(pos):<3d}              {n:3d}/{len(neg):<3d}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tables", required=True, help="glob for *.offtarget_analysis.tsv")
    ap.add_argument("--queue-all", required=True, help="review_queue_all.tsv")
    ap.add_argument("--truth-wgs", help="cart_wgs_merged.xlsx")
    ap.add_argument("--figdir")
    a = ap.parse_args()

    paths = sorted(glob.glob(a.tables))
    df = nm.load_tables(paths)
    ev = df[df.indel_reads > 0]
    print(f"tables {len(paths)}   site-rows {len(df)}   rows with indel evidence {len(ev)}")

    # ---------------------------------------------------------------------------------------
    # 1. Background: is it Uniform?
    # ---------------------------------------------------------------------------------------
    print("\n" + "=" * 94)
    print("BACKGROUND — off-target rows with indel evidence")
    print("=" * 94)
    bg = hist(ev[ev.is_target == 0].min_cut_distance)
    pct = bg / bg.sum() * 100
    print(f"n = {bg.sum()}     bins 0..{WINDOW-1}")
    print("counts :", bg.tolist())
    print("percent:", np.round(pct, 2).tolist(), f"   uniform expectation {100/WINDOW:.2f}")
    chi = chisquare(bg)
    print(f"chi-square vs Uniform: stat {chi.statistic:.1f}, p = {chi.pvalue:.4f}")
    print(f"spread: min {pct.min():.2f}%  max {pct.max():.2f}%  "
          f"(largest deviation from uniform {np.abs(pct - 100/WINDOW).max():.2f} points)")
    if chi.pvalue < 0.05:
        print("NOTE: strict uniformity is rejected. At n = %d even a fraction-of-a-point wobble is"
              % bg.sum())
        print("      detectable, and every bin is within ~1.1 points of flat, so Uniform remains a")
        print("      good working null -- but it is an approximation, not an exact property.")

    # ---------------------------------------------------------------------------------------
    # 2. Signal
    # ---------------------------------------------------------------------------------------
    print("\n" + "=" * 94)
    print("SIGNAL — on-target rows with indel evidence (real Cas9 cuts)")
    print("=" * 94)
    on = ev[ev.is_target == 1]
    sg = hist(on.min_cut_distance)
    print(f"n = {sg.sum()}")
    print("counts :", sg.tolist())
    print("percent:", np.round(sg / sg.sum() * 100, 2).tolist())
    print(f"mass at d <= 1: {sg[:2].sum()}/{sg.sum()} = {sg[:2].sum()/sg.sum()*100:.1f}%")
    dvals = pd.to_numeric(on.min_cut_distance, errors="coerce").dropna().astype(int)
    inw = dvals[dvals < WINDOW]
    sigma = fit_halfnormal(inw)
    print(f"discretised half-normal MLE: sigma = {sigma:.3f} bp")
    print("half-normal pmf:", np.round(halfnormal_pmf(sigma), 4).tolist())
    # sqrt(mean(d^2)) is driven by the largest observation, so with n this small a single stray
    # event moves it a long way. Say so rather than presenting the fit as if it were stable.
    sig_trim = fit_halfnormal(inw[inw <= 2])
    print(f"same fit excluding events beyond 2 bp ({int((inw>2).sum())} of {len(inw)}): "
          f"sigma = {sig_trim:.3f} bp")
    print(f"The MLE is not robust here: it puts {halfnormal_pmf(sigma)[0]:.2f} at d=0 against an")
    print(f"empirical {sg[0]/sg.sum():.2f}, because sqrt(mean(d^2)) is dominated by the tail. The")
    print("LR below therefore uses the smoothed EMPIRICAL pmf; the half-normal is reported only as")
    print("a transferable one-parameter summary, and should be refitted per cohort if used.")

    # ---------------------------------------------------------------------------------------
    # 3. The likelihood ratio
    # ---------------------------------------------------------------------------------------
    print("\n" + "=" * 94)
    print("LIKELIHOOD RATIO  LR(d) = P(d|signal) / P(d|background)")
    print("=" * 94)
    p_sig = smoothed_pmf(sg)
    p_bg = np.full(WINDOW, 1.0 / WINDOW)
    lr = p_sig / p_bg
    print("  d   P(d|signal)   P(d|bg)     LR      10log10 LR (evidence, dB)")
    for d in range(WINDOW):
        print(f" {d:2d}      {p_sig[d]:.4f}      {p_bg[d]:.4f}   {lr[d]:7.3f}   "
              f"{10*np.log10(lr[d]):+7.2f}")
    print("\nThe last column is on the same decibel scale as AQ, which is the point: a distance of")
    print(f"0 contributes {10*np.log10(lr[0]):+.1f} dB of evidence and a distance of 10 contributes")
    print(f"{10*np.log10(lr[-1]):+.1f} dB, instead of both being 'inside the window, therefore fine'.")

    # ---------------------------------------------------------------------------------------
    # 4. Validation against the curated label
    # ---------------------------------------------------------------------------------------
    if a.truth_wgs:
        lab = nm.wgs_curated_label(a.truth_wgs)
        q = pd.read_csv(a.queue_all, sep="\t")
        q["guide"] = q.sample_name.map(nm.guide_of)
        q["chrom"] = q.chrom.astype(str)
        m = q.merge(lab, on=["guide", "chrom", "start"], how="inner")
        pos = pd.to_numeric(m[m.curated_label == 1].cut_dist_min, errors="coerce")
        neg = pd.to_numeric(m[m.curated_label == 0].cut_dist_min, errors="coerce")
        print("\n" + "=" * 94)
        print("VALIDATION — curated WGS label, gated rows")
        print("=" * 94)
        print(f"joined {len(m)} rows: {len(pos)} confirmed / {len(neg)} human-rejected")
        sweep(pos, neg)
        safe = [t for t in range(WINDOW) if int((pos <= t).sum()) == len(pos)]
        if safe:
            t0 = min(safe)
            print(f"\ntightest threshold holding full recall: d <= {t0}")
            print(f"  at d <= {t0}: keeps {int((neg<=t0).sum())}/{len(neg)} rejected rows, "
                  f"vs {int((neg<=10).sum())}/{len(neg)} at the shipped d <= 10")

        # What it is worth AFTER the other rules have run -- the only number that decides anything.
        kept = m[m.why_dropped.fillna("") == ""]
        allq = q[q.why_dropped.fillna("") == ""]
        kp = pd.to_numeric(kept[kept.curated_label == 1].cut_dist_min, errors="coerce")
        kn = pd.to_numeric(kept[kept.curated_label == 0].cut_dist_min, errors="coerce")
        aq = pd.to_numeric(allq.cut_dist_min, errors="coerce")
        print("\n-- marginal effect on the FINAL queue, after every other rule has run --")
        print(f"queue {len(allq)} rows; labelled subset {len(kp)} confirmed / {len(kn)} rejected")
        for t in (1, 2, 3, 5, 10):
            print(f"  d<={t:2d}: queue {int((aq<=t).sum()):3d}/{len(allq)}   "
                  f"confirmed {int((kp<=t).sum())}/{len(kp)}   rejected {int((kn<=t).sum())}/{len(kn)}")
        print("\nIf those rows barely move, the spatial rule is REDUNDANT with the rules ahead of")
        print("it on this cohort, and its value is as a score for ranking and for the single-sample")
        print("case where no matched control exists -- not as a tighter gate. Report it that way.")

    # ---------------------------------------------------------------------------------------
    # 5. Tail, from the uncapped anchor distance -- a DIFFERENT quantity
    # ---------------------------------------------------------------------------------------
    print("\n" + "=" * 94)
    print("TAIL — from indel_info field 7 (anchor distance, uncapped). NOT min_cut_distance.")
    print("=" * 94)
    tail = []
    for info in ev.indel_info.dropna():
        for e in str(info).split(";"):
            f = e.split("|")
            if len(f) >= 11:
                try:
                    tail.append(abs(int(f[7])))
                except ValueError:
                    pass
    tail = np.array(tail)
    if len(tail):
        bc = np.bincount(tail)
        print(f"events {len(tail)}   median {np.median(tail):.0f}   max {tail.max()}")
        print("counts 0..:", bc.tolist())
        print(f"fraction beyond 10 bp: {(tail>10).mean()*100:.2f}%")
        print("The anchor distance is always >= min_cut_distance, so this over-states the tail of")
        print("the quantity the filter uses. It shows the tail EXISTS and is thin; it does not")
        print("give its shape for min_cut_distance. Only a re-run with a larger -d could do that.")

    # ---------------------------------------------------------------------------------------
    if a.figdir:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        os.makedirs(a.figdir, exist_ok=True)
        fig, ax = plt.subplots(1, 3, figsize=(14, 4))
        x = np.arange(WINDOW)
        ax[0].bar(x, bg / bg.sum(), color="#777")
        ax[0].axhline(1 / WINDOW, ls="--", color="crimson", label="Uniform 1/11")
        ax[0].set_title(f"background: off-target (n={bg.sum()})")
        ax[0].set_xlabel("distance from predicted cut (bp)"); ax[0].legend()
        ax[1].bar(x, sg / sg.sum(), color="#2a6")
        ax[1].plot(x, halfnormal_pmf(sigma), "o--", color="k",
                   label=f"half-normal $\\sigma$={sigma:.2f}")
        ax[1].set_title(f"signal: on-target (n={sg.sum()})")
        ax[1].set_xlabel("distance from predicted cut (bp)"); ax[1].legend()
        ax[2].bar(x, 10 * np.log10(lr), color="#36c")
        ax[2].axhline(0, color="k", lw=0.8)
        ax[2].set_title("evidence  $10\\log_{10}$ LR(d)")
        ax[2].set_xlabel("distance from predicted cut (bp)"); ax[2].set_ylabel("dB")
        fig.tight_layout()
        p = os.path.join(a.figdir, "cut_distance_model.png")
        fig.savefig(p, dpi=150)
        print(f"\nwrote {p}")


if __name__ == "__main__":
    main()
