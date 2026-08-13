#!/usr/bin/env python3
"""
noise_model.py — score off-target calls against a probabilistic noise baseline.

Why this exists
---------------
`review_filter.py` treats every noise source as a boolean: rule 4 asks "is this site in the panel
of normals?", rule 6 asks "does this site touch the DRAGEN systematic-noise panel?". A boolean is
the wrong shape for the question. A site can carry 1% background in the controls and 40% signal in
the treated sample; existence-based filtering throws that edit away, and the measured cost is real
(a bare panel-interval hit flags 33.1% of artifacts but also 1.2% of genuine on-target edits).

This script replaces existence with expectation. For a call with `k` indel reads out of `n`, it
asks how surprising `k` is under a per-locus background model, and reports the answer Phred-scaled:

    AQ = -10 log10 P(X >= k | n, background at this locus)

High AQ = the observation cannot be explained by background = keep. Low AQ = indistinguishable from
this locus's known noise = drop.

Why beta-binomial and not binomial
----------------------------------
The noise is overdispersed, and demonstrably so: 61.3% of loci in the DRAGEN IDPF panel are flagged
in exactly ONE panel sample. A point-`p` binomial assumes every sample sees the same error rate, so
at a locus where 1 of 46 samples showed 17% and the other 45 showed nothing, it uses p = 0.0037 and
declares a 5% observation overwhelmingly significant. The beta-binomial carries the spread as well
as the mean, which is the only reason the tail probability means anything -- and the tail is the
entire output.

Baselines (--baseline)
----------------------
  matched   global prior + THIS sample's own matched control at the locus.  Replaces rule 1:
            a germline het sits near VAF 0.5 in the matched control, so a 0.5 observation in the
            treated sample is unsurprising and scores low without a separate germline rule.
  loo       global prior + the OTHER samples' controls at the locus (leave-one-out). Replaces
            rule 4: this is cohort recurrence expressed as a rate rather than a blacklist.
  both      both sets pooled.
  panel     the external DRAGEN panel, no controls at all -- the single-sample case. Combined
            with --panel-p below.

Empirical Bayes, which is also how the "no data here" case is handled
---------------------------------------------------------------------
A Beta(a0, b0) prior is fitted once by method of moments over every control observation in the
input, then updated per locus:  a = a0 + (control alt reads), b = b0 + (control ref reads). The
predictive distribution for a new sample is BetaBinom(n, a, b).

This is what makes the floor principled rather than arbitrary. ~97% of genomic positions are absent
from the DRAGEN panel (it covers 3.16% of the genome), and a missing locus is NOT p=0 -- that would
make every observation infinitely significant. Here an unobserved locus simply keeps the global
prior, which is the honest statement that we know nothing beyond the assay-wide error rate.

--panel-p, and why 'mean' is wrong
-----------------------------------
The panel's MEAN column is diluted by the whole panel: verified MAX/MEAN = 45.99 at NR=1, exactly
the 46-sample panel size, so MEAN = (sum of VAFs) / N_panel and includes the ~45 samples that showed
nothing. Using it as the background rate understates the rate an AFFECTED sample sees by up to N/NR.
  mean       p = MEAN                 (anti-conservative; provided for comparison)
  max        p = MAX                  (an order statistic over N draws; conservative)
  corrected  p = MEAN * N / NR        (mean rate among the samples that actually showed noise)
N is read from the panel's own `##PON SAMPLES:` header.

usage:
  noise_model.py IN.tsv [IN2.tsv ...] -o scored.tsv --baseline loo
  noise_model.py IN.tsv ... -o scored.tsv --baseline panel \
                 --snv-noise IDPF_WGS_hg38_v.2.0.0_systematic_noise.snv.bed.gz --panel-p corrected
"""
import argparse
import collections
import glob
import gzip
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import betabinom, binom

SNV_NOISE_SLOP = 2
# Floor on the prior mean. Without it a cohort with pristine controls yields a0/(a0+b0) ~ 0 and
# every call becomes infinitely significant -- the same failure as p=0 for an absent locus.
MIN_PRIOR_MEAN = 1e-5


def sample_name_of(path):
    return os.path.basename(path).split(".")[0]


def load_tables(paths):
    frames = []
    for p in paths:
        df = pd.read_csv(p, sep="\t")
        df["sample_name"] = sample_name_of(p)
        frames.append(df)
    df = pd.concat(frames, ignore_index=True)
    for c in ("total_reads", "indel_reads", "control_reads", "control_indel_reads"):
        df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0)
    return df


def fit_global_prior(alt, depth):
    """Method-of-moments Beta prior over per-observation control VAFs.

    Only observations with depth > 0 carry information. The variance is taken across observations,
    so it absorbs both sampling noise and true locus-to-locus rate variation -- which is what we
    want the prior to express.
    """
    ok = depth > 0
    if ok.sum() < 2:
        return 0.5, 500.0
    vaf = alt[ok] / depth[ok]
    m = float(np.mean(vaf))
    v = float(np.var(vaf, ddof=1))
    m = max(m, MIN_PRIOR_MEAN)
    if v <= 0 or v >= m * (1 - m):
        # degenerate: fall back to a weakly-informative prior centred on the observed mean
        return m * 100.0, (1 - m) * 100.0
    k = m * (1 - m) / v - 1
    return max(m * k, 1e-6), max((1 - m) * k, 1e-6)


def locus_control_support(df):
    """(chrom, end) -> [(sample, control_alt, control_depth), ...]"""
    d = collections.defaultdict(list)
    for s, c, e, a, n in zip(df.sample_name, df.chrom, df.end,
                             df.control_indel_reads, df.control_reads):
        d[(c, int(e))].append((s, float(a), float(n)))
    return d


def load_panel(path, positions, slop=SNV_NOISE_SLOP):
    """DRAGEN SNV panel -> {(chrom,pos): (mean, max, alleles, nr)} plus the panel sample count.

    Streams the file; only queried positions are retained. The `##PON SAMPLES:` header gives N,
    which --panel-p corrected needs to undo the MEAN dilution.
    """
    want = collections.defaultdict(set)
    for c, e in positions:
        for d in range(-slop, slop + 1):
            want[c].add(int(e) + d)
    hits, n_panel = {}, None
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as f:
        for line in f:
            if line[0] == "#":
                if line.startswith("##PON SAMPLES:"):
                    n_panel = len([x for x in line.split(":", 1)[1].split(",") if x.strip()])
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 7:
                continue
            s = want.get(p[0])
            if s is None:
                continue
            try:
                pos = int(p[2])
                if pos not in s:
                    continue
                rec = (float(p[3]), float(p[4]), p[5], int(p[6]))
            except ValueError:
                continue
            prev = hits.get((p[0], pos))
            if prev is None or rec[3] > prev[3]:
                hits[(p[0], pos)] = rec
    return hits, (n_panel or 46)


def panel_rate(hits, chrom, end, n_panel, mode, floor_p, slop=SNV_NOISE_SLOP):
    """Background rate at one locus from the panel, or the floor when the panel says nothing.

    Only indel-capable records count: a locus that is noisy for substitutions says nothing about
    an indel call at the same coordinate.
    """
    best = None
    for d in range(-slop, slop + 1):
        rec = hits.get((chrom, int(end) + d))
        if rec and ("D" in rec[2] or "I" in rec[2]):
            if best is None or rec[3] > best[3]:
                best = rec
    if best is None:
        return floor_p, 0
    mean, mx, _alleles, nr = best
    if mode == "mean":
        p = mean
    elif mode == "max":
        p = mx
    else:
        p = mean * n_panel / max(nr, 1)
    return min(max(p, floor_p), 0.999), nr


def control_posterior(df, q, baseline, a0, b0):
    """Per-row Beta posterior over the background rate, from control observations.

    Returns (alpha, beta, n_observations, pooled_control_depth). The depth is returned because
    it bounds how much the posterior is entitled to claim -- see apply_depth_floor.
    """
    support = locus_control_support(df)
    alpha = np.full(len(q), a0, dtype=float)
    beta = np.full(len(q), b0, dtype=float)
    nsup = np.zeros(len(q), dtype=int)
    dep_out = np.zeros(len(q), dtype=float)
    for i, (s, c, e) in enumerate(zip(q.sample_name, q.chrom, q.end)):
        alt = dep = 0.0
        cnt = 0
        for samp, ca, cn in support.get((c, int(e)), ()):
            if cn <= 0:
                continue
            if baseline == "matched" and samp != s:
                continue
            if baseline == "loo" and samp == s:
                continue
            alt += ca; dep += cn; cnt += 1
        alpha[i] = a0 + alt
        beta[i] = b0 + max(dep - alt, 0.0)
        nsup[i] = cnt
        dep_out[i] = dep
    return alpha, beta, nsup, dep_out


def apply_depth_floor(alpha, beta, control_depth):
    """Stop the posterior claiming more resolution than the control actually has.

    A control with zero alt reads at depth d does NOT show the background is ~0. It shows the
    background is below roughly 1/d, and nothing more. Left alone the posterior says otherwise:
    with alpha = a0 and beta = b0 + d, the mean collapses to a0/(a0+d) -- on this cohort 3.5e-5
    against a 170x control, about 70x below the fitted prior and ~100x below what 170 reads can
    support. That is the classic clean-control trap, and it inflates AQ for every call at a locus
    the control simply never sampled deeply enough to speak about.

    So floor the posterior MEAN at 1/control_depth, holding the concentration (alpha+beta) fixed
    so only the location moves and the strength of belief is preserved. Rows with no control depth
    are left on the global prior -- there is no depth to derive a floor from.

    Today the VAF>=0.005 gate hides this (only 2 gated rows sit below a 170x control's 0.59%
    resolution). It stops being hidden the moment anyone lowers the gate for a high-sensitivity run.
    """
    alpha = np.asarray(alpha, dtype=float).copy()
    beta = np.asarray(beta, dtype=float).copy()
    depth = np.asarray(control_depth, dtype=float)
    conc = alpha + beta
    mean = alpha / conc
    has_depth = depth > 0
    floor = np.where(has_depth, 1.0 / np.maximum(depth, 1.0), 0.0)
    lift = has_depth & (mean < floor)
    alpha[lift] = floor[lift] * conc[lift]
    beta[lift] = conc[lift] - alpha[lift]
    return alpha, beta, int(lift.sum())


def aq_from_sf(sf):
    """Phred-scale a survival probability, capped so log10(0) does not become inf."""
    sf = np.clip(np.asarray(sf, dtype=float), 1e-300, 1.0)
    return -10.0 * np.log10(sf)


# ---------------------------------------------------------------------------------------------
# Truth. Two sources, deliberately kept separate because they answer different questions.
# ---------------------------------------------------------------------------------------------
GUIDE_ALIAS = {"CTLA41": "CTLA4"}
CHROM_REPAIR = {"c": "chr5"}          # one corrupted cell in CART_NS0011-CREBRF (chr5q35.1)
CONFIRMED = {"1", "1.0", "1?"}
# The WGS review adjudicated every row inside this stratum by eye, so a blank manual_review
# there means REJECTED. Outside it, blank means unreviewed and the row carries no label.
WGS_STRATUM_IF, WGS_STRATUM_READS = 0.05, 10


def guide_of(name):
    """Pipeline sample or review-sheet name -> guide.

    `CART_NS0011-ABTB1` -> ABTB1, `CART_NS0027-CTLA41_2` -> CTLA4, `ABTB1-KO-DNA` -> ABTB1,
    `ATF7IP-KO-CART-DNA` -> ATF7IP.
    """
    s = str(name)
    import re
    if re.match(r"^(CART_)?NS\d+-", s):
        g = s.split("-", 1)[1]
    else:
        g = s.split("-", 1)[0]
    g = re.sub(r"_\d+$", "", g)
    return GUIDE_ALIAS.get(g, g)


def wgs_curated_label(xlsx, sheet="gold_wgs"):
    """(guide, chrom, start) -> 1 confirmed / 0 human-rejected, inside the reviewed stratum only."""
    g = pd.read_excel(xlsx, sheet_name=sheet)
    iff = pd.to_numeric(g.indel_fraction, errors="coerce")
    ir = pd.to_numeric(g.indel_reads, errors="coerce")
    strat = (iff >= WGS_STRATUM_IF) & (ir >= WGS_STRATUM_READS)
    conf = g["manual_review"].astype(str).str.strip().isin(CONFIRMED)
    r = g[strat].copy()
    r["curated_label"] = conf[strat].astype(int)
    r["guide"] = r["sample_name"].map(guide_of)
    r["chrom"] = r["chrom"].astype(str).replace(CHROM_REPAIR)
    return r.groupby(["guide", "chrom", "start"], as_index=False)["curated_label"].max()


def ecs_confirmed(csv):
    """(guide, chrom, start) for edits a human confirmed in the ECS assay.

    Positives only -- NaN in that table means NOT REVIEWED, not rejected -- so this supports
    recall and nothing else. Its value is reach: ECS runs at ~1,960x, so it confirms edits well
    below the 5% floor of the WGS review.
    """
    g = pd.read_csv(csv, low_memory=False)
    g = g[g["manual_review"].astype(str).str.strip().isin(CONFIRMED)].copy()
    g["guide"] = g["sample_name"].map(guide_of)
    g["chrom"] = g["chrom"].astype(str).replace(CHROM_REPAIR)
    g["ecs_vaf"] = pd.to_numeric(g["indel_fraction"], errors="coerce")
    return g[["guide", "chrom", "start", "ecs_vaf"]].drop_duplicates(["guide", "chrom", "start"])


def evaluate(q, truth_wgs, truth_ecs):
    q = q.copy()
    q["guide"] = q.sample_name.map(guide_of)
    q["chrom"] = q.chrom.astype(str)

    if truth_wgs:
        lab = wgs_curated_label(truth_wgs)
        m = q.merge(lab, on=["guide", "chrom", "start"], how="inner")
        pos, neg = m[m.curated_label == 1], m[m.curated_label == 0]
        print(f"\n== curated WGS label (two-class, VAF>=5% stratum) ==")
        print(f"joined {len(m)} scored rows: {len(pos)} confirmed / {len(neg)} human-rejected")
        if len(pos) and len(neg):
            from sklearn.metrics import roc_auc_score, average_precision_score
            y, s = m.curated_label.values, m.AQ.values
            print(f"AUC {roc_auc_score(y, s):.3f}   AP {average_precision_score(y, s):.3f}")
            print(f"AQ  confirmed: median {pos.AQ.median():.1f}  min {pos.AQ.min():.1f}")
            print(f"AQ  rejected : median {neg.AQ.median():.1f}  max {neg.AQ.max():.1f}")
            for thr in (10, 20, 30, 60):
                tp = int((pos.AQ >= thr).sum()); fp = int((neg.AQ >= thr).sum())
                prec = tp / (tp + fp) if tp + fp else float("nan")
                print(f"  AQ>={thr:<3d} recall {tp}/{len(pos)} = {tp/len(pos):.3f}"
                      f"   precision {prec:.3f}   (FP {fp})")

    if truth_ecs:
        e = ecs_confirmed(truth_ecs)
        m = q.merge(e, on=["guide", "chrom", "start"], how="inner")
        sub = m[m.ecs_vaf < 0.05]
        print(f"\n== ECS-confirmed edits (recall only; positives-only truth) ==")
        print(f"joined {len(m)} confirmed edits, {len(sub)} of them below 5% ECS VAF")
        for thr in (10, 20, 30, 60):
            print(f"  AQ>={thr:<3d} all {int((m.AQ>=thr).sum())}/{len(m)}"
                  + (f"   sub-5% {int((sub.AQ>=thr).sum())}/{len(sub)}" if len(sub) else ""))
        if len(sub):
            print("  sub-5% edits (ECS VAF -> WGS AQ): "
                  + ", ".join(f"{v:.3f}->{aq:.1f}" for v, aq in zip(sub.ecs_vaf, sub.AQ)))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("inputs", nargs="+", help="*.offtarget_analysis.tsv")
    ap.add_argument("-o", "--out", required=True)
    ap.add_argument("--baseline", default="loo", choices=["matched", "loo", "both", "panel"])
    ap.add_argument("--snv-noise", metavar="BED", help="DRAGEN SNV panel (--baseline panel)")
    ap.add_argument("--panel-p", default="corrected", choices=["mean", "max", "corrected"])
    ap.add_argument("--floor-p", type=float, default=1e-3,
                    help="background rate for loci the panel does not mention (default 1e-3). "
                         "~97%% of positions are absent from the panel, so this governs most sites")
    ap.add_argument("--min-reads", type=int, default=2)
    ap.add_argument("--min-vaf", type=float, default=0.005)
    ap.add_argument("--no-depth-floor", dest="depth_floor", action="store_false",
                    help="disable the 1/control_depth floor on the posterior background. The "
                         "floor is on by default; turning it off restores the clean-control trap "
                         "and is provided only to measure its effect")
    ap.add_argument("--truth-wgs", metavar="XLSX",
                    help="cart_wgs_merged.xlsx -- two-class curated label, VAF>=5% stratum only")
    ap.add_argument("--truth-ecs", metavar="CSV",
                    help="cart_ecs_merged.csv.gz -- confirmed edits incl. sub-5%% (recall only)")
    a = ap.parse_args()

    paths = []
    for p in a.inputs:
        paths.extend(sorted(glob.glob(p)) if any(ch in p for ch in "*?[") else [p])
    df = load_tables(paths)
    n_in = len(df)

    vaf = df.indel_reads / df.total_reads.replace(0, np.nan)
    q = df[(df.indel_reads >= a.min_reads) & (vaf.fillna(0) >= a.min_vaf)].copy().reset_index(drop=True)
    print(f"input rows            : {n_in}")
    print(f"cleared the gate      : {len(q)}   (reads>={a.min_reads}, VAF>={a.min_vaf})")

    a0, b0 = fit_global_prior(df.control_indel_reads.values, df.control_reads.values)
    print(f"global prior          : Beta(a={a0:.4g}, b={b0:.4g})  mean={a0/(a0+b0):.6f}")

    k = q.indel_reads.values.astype(int)
    n = q.total_reads.values.astype(int)

    if a.baseline == "panel":
        if not a.snv_noise:
            sys.exit("ERROR: --baseline panel requires --snv-noise")
        hits, n_panel = load_panel(a.snv_noise, set(zip(q.chrom, q.end)))
        print(f"panel                 : {len(hits)} records at queried positions, N={n_panel}")
        ps, nrs = [], []
        for c, e in zip(q.chrom, q.end):
            p, nr = panel_rate(hits, c, e, n_panel, a.panel_p, a.floor_p)
            ps.append(p); nrs.append(nr)
        ps = np.asarray(ps)
        q["bg_rate"] = ps
        q["panel_nr"] = nrs
        q["from_panel"] = np.asarray(nrs) > 0
        # a point rate is all the panel gives, so this arm is an honest binomial
        q["AQ"] = aq_from_sf(binom.sf(k - 1, n, ps))
        print(f"loci with panel support: {int(q.from_panel.sum())} / {len(q)} "
              f"({q.from_panel.mean()*100:.1f}%) -- the rest use --floor-p {a.floor_p}")
    else:
        alpha, beta, nsup, dep = control_posterior(df, q, a.baseline, a0, b0)
        if a.depth_floor:
            alpha, beta, n_floored = apply_depth_floor(alpha, beta, dep)
            print(f"depth floor           : raised background on {n_floored} rows "
                  f"({n_floored/len(q)*100:.1f}%) to 1/control_depth")
        q["bg_rate"] = alpha / (alpha + beta)
        q["n_control_obs"] = nsup
        q["control_depth"] = dep
        q["AQ"] = aq_from_sf(betabinom.sf(k - 1, n, alpha, beta))
        print(f"control observations  : median {int(np.median(nsup))} per locus "
              f"({int((nsup == 0).sum())} loci with none -> global prior)")

    q["AQ"] = q["AQ"].round(2)
    q["bg_rate"] = q["bg_rate"].round(6)
    q.to_csv(a.out, sep="\t", index=False)

    for thr in (10, 20, 30, 60):
        n_pass = int((q.AQ >= thr).sum())
        print(f"  AQ >= {thr:<3d} : {n_pass:5d} rows kept "
              f"(chance FPs at {len(q)} tests: {len(q) * 10 ** (-thr / 10):.1f})")
    print(f"scored {len(q)} rows -> {a.out}")

    if a.truth_wgs or a.truth_ecs:
        evaluate(q, a.truth_wgs, a.truth_ecs)


if __name__ == "__main__":
    main()
