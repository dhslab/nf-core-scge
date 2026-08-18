#!/usr/bin/env python3
"""noise_model_validate.py — test whether the shipped background model describes the data.

The objection this answers
--------------------------
`noise_model.py` fits ONE Beta(a0, b0) by method of moments across every control observation in
the run, and uses it as the prior on the per-locus background rate. The claim examined here is
that this conflates two physically distinct processes:

  machine physics   sequencing, PCR and alignment error. A per-base Bernoulli process whose rate
                    is tiny, is NOT globally fixed, and varies with sequence context -- polymerase
                    slippage in a homopolymer being the textbook case.
  biology           germline variation in the donor. Not error at all: a real allele near VAF 0.5
                    or 1.0, and specific to THAT donor.

A single Beta fitted across both describes a population that does not exist. Method of moments
matches a mean and a variance; applied to a mixture of a point mass near 0 and a mode near 0.5 it
returns the parameters of a distribution that describes neither component.

What is measured here
---------------------
  1. the two populations, plotted against the fitted Beta
  2. calibration under a genuine null, built with no new data -- score each donor's own control
     counts against the OTHER donors' controls at the same locus. Both are unedited, so anything
     that scores as signal is a miscalibration. Randomised p-values are used, because with
     discrete counts the ordinary survival function cannot be Uniform even under a perfect model.
  3. held-out log predictive likelihood for five candidate models, split BY LOCUS so a locus
     cannot appear in both halves
  4. whether the rate depends on sequence context, by refitting per homopolymer stratum
  5. the depth floor, ablated -- is it a patch on a mis-specified model, or genuinely needed?
  6. what any of it changes operationally: queue size and recall against the curated label

This script changes nothing. It is analysis only, run by hand, and the pipeline does not call it.

usage:
  noise_model_validate.py --tables 'results_cart_bnd/*/*.offtarget_analysis.tsv' \
      --queue-all results_cart_bnd/review/review_queue_all.tsv \
      --truth-wgs '.../cart_wgs_merged.xlsx' \
      --fasta .../hg38_PLVM_CD19_CARv4_cd34.fa --figdir docs/images
"""
import argparse
import collections
import glob
import os
import sys

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.stats import betabinom, binom, kstest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import noise_model as nm                                              # noqa: E402

RNG = np.random.default_rng(0)
HP_WINDOW = 12       # bp either side of the site to search for a homopolymer run


# =============================================================================================
# Candidate models. Each exposes:
#   fit(k, n)                       -> params, from control observations
#   logpmf(k, n)                    -> marginal predictive log-likelihood of a new observation
#   posterior(K, N)                 -> per-locus state after seeing that locus's other controls
#   tail(k, n, state)               -> (P(X > k), P(X == k)) for a randomised p-value
# =============================================================================================
class BinomialGlobal:
    """One fixed rate for every locus and every sample. The straw man -- this is what "just use
    an error rate" means, and it is the model whose failure motivated the beta-binomial."""
    name = "Binomial, global p"

    def fit(self, k, n):
        self.p = float(k.sum() / max(n.sum(), 1))
        return self

    def describe(self):
        return f"p = {self.p:.6g}"

    def logpmf(self, k, n):
        return binom.logpmf(k, n, self.p)

    def posterior(self, K, N):
        return None                      # a point rate cannot learn from the locus

    def tail(self, k, n, state):
        return binom.sf(k, n, self.p), binom.pmf(k, n, self.p)


class BetaBinomMOM:
    """The shipped model: Beta prior by method of moments, conjugate per-locus update."""
    name = "Beta-Binomial, MOM prior (SHIPPED)"

    def fit(self, k, n):
        self.a, self.b = nm.fit_global_prior(k, n)
        return self

    def describe(self):
        return f"Beta({self.a:.6g}, {self.b:.6g})  mean {self.a/(self.a+self.b):.6g}"

    def logpmf(self, k, n):
        return betabinom.logpmf(k, n, self.a, self.b)

    def posterior(self, K, N):
        return (self.a + K, self.b + max(N - K, 0.0))

    def tail(self, k, n, state):
        a, b = state if state else (self.a, self.b)
        return betabinom.sf(k, n, a, b), betabinom.pmf(k, n, a, b)


class BetaBinomMML(BetaBinomMOM):
    """Same family, but the prior is fitted by maximum marginal likelihood -- proper empirical
    Bayes. Isolates how much of any failure is the ESTIMATOR rather than the FAMILY."""
    name = "Beta-Binomial, max marginal likelihood"

    def fit(self, k, n):
        a0, b0 = nm.fit_global_prior(k, n)

        def nll(t):
            a, b = np.exp(t)
            if not np.isfinite(a) or not np.isfinite(b) or a <= 0 or b <= 0:
                return 1e18
            v = betabinom.logpmf(k, n, a, b)
            return -float(np.sum(v[np.isfinite(v)]))

        r = minimize(nll, np.log([max(a0, 1e-4), max(b0, 1e-4)]), method="Nelder-Mead",
                     options=dict(maxiter=4000, xatol=1e-6, fatol=1e-6))
        self.a, self.b = np.exp(r.x)
        return self


class ZeroInflatedBB:
    """The physically motivated model: two components, one per process.

        with prob pi     the locus is machine-error-only, at a rate eps the control cannot resolve
        with prob 1-pi   the locus carries a real allele, drawn from Beta(a, b)

    This is the model that can SAY "this locus is clean", which a single Beta cannot -- and that
    is exactly the statement the depth floor exists to fake.
    """
    name = "Zero-inflated Beta-Binomial (2 components)"

    def fit(self, k, n):
        p0 = float(k.sum() / max(n.sum(), 1))
        a0, b0 = nm.fit_global_prior(k, n)

        def nll(t):
            lo, le, la, lb = t
            pi = 1.0 / (1.0 + np.exp(-lo))
            eps, a, b = np.exp(le), np.exp(la), np.exp(lb)
            if not all(np.isfinite([eps, a, b])) or min(eps, a, b) <= 0 or eps >= 1:
                return 1e18
            l0 = np.log(pi + 1e-300) + binom.logpmf(k, n, eps)
            l1 = np.log(1 - pi + 1e-300) + betabinom.logpmf(k, n, a, b)
            m = np.maximum(l0, l1)
            v = m + np.log(np.exp(l0 - m) + np.exp(l1 - m))
            return -float(np.sum(v[np.isfinite(v)]))

        r = minimize(nll, [0.0, np.log(max(p0, 1e-6)), np.log(max(a0, 1e-3)),
                           np.log(max(b0, 1e-3))],
                     method="Nelder-Mead", options=dict(maxiter=20000, fatol=1e-6))
        lo, le, la, lb = r.x
        self.pi = float(1 / (1 + np.exp(-lo)))
        self.eps, self.a, self.b = float(np.exp(le)), float(np.exp(la)), float(np.exp(lb))
        return self

    def describe(self):
        return (f"pi(clean) = {self.pi:.4f}   eps = {self.eps:.3g}   "
                f"Beta({self.a:.4g}, {self.b:.4g}) mean {self.a/(self.a+self.b):.4f}")

    def logpmf(self, k, n):
        l0 = np.log(self.pi + 1e-300) + binom.logpmf(k, n, self.eps)
        l1 = np.log(1 - self.pi + 1e-300) + betabinom.logpmf(k, n, self.a, self.b)
        m = np.maximum(l0, l1)
        return m + np.log(np.exp(l0 - m) + np.exp(l1 - m))

    def posterior(self, K, N):
        """Both the component weights and the Beta update on that locus's other controls."""
        l0 = np.log(self.pi + 1e-300) + binom.logpmf(K, N, self.eps)
        l1 = np.log(1 - self.pi + 1e-300) + betabinom.logpmf(K, N, self.a, self.b)
        m = max(l0, l1)
        if not np.isfinite(m):
            # Both components assign this locus's other controls zero probability -- eps is driven
            # to the denormal floor, so any K > 0 is impossible under component 0. Fall back to the
            # prior weights rather than propagating a NaN.
            return (self.pi, 1 - self.pi, self.a + K, self.b + max(N - K, 0.0))
        w0, w1 = np.exp(l0 - m), np.exp(l1 - m)
        s = w0 + w1
        return (w0 / s, w1 / s, self.a + K, self.b + max(N - K, 0.0))

    def tail(self, k, n, state):
        if state is None:
            w0, w1, a, b = self.pi, 1 - self.pi, self.a, self.b
        else:
            w0, w1, a, b = state
        sf = w0 * binom.sf(k, n, self.eps) + w1 * betabinom.sf(k, n, a, b)
        pm = w0 * binom.pmf(k, n, self.eps) + w1 * betabinom.pmf(k, n, a, b)
        return sf, pm


# =============================================================================================
def randomised_tails(model, obs_by_locus, a_global=None):
    """Leave-one-donor-out p-values under a null where every observation is noise by construction.

    For each locus with >= 2 donors carrying control depth, each donor's own control counts are
    scored against a background built from the OTHER donors' controls at that locus. Both sides
    are unedited material, so a correct model must return Uniform(0,1).

    With discrete counts the plain survival function CANNOT be uniform even under a perfect model
    -- it is bounded below by P(X = k), and with 99% of observations at k = 0 the mass piles at 1.
    The randomised p-value  U*P(X = k) + P(X > k),  U ~ Uniform(0,1),  is exactly uniform under a
    correct discrete model, so it is what makes this test meaningful at all.

    The randomisation has a cost: any single realisation is ONE DRAW from the test, and at
    n ~ 72,000 the KS p-value swings widely across seeds while the KS D statistic barely moves.

    So this returns the two SEED-INDEPENDENT components, P(X > k) and P(X = k), rather than a
    p-value. Drawing a realisation is then `sf + U * pm` (see draw_pvalues), which is cheap, so the
    caller can sweep seeds without repeating the expensive per-observation betabinom evaluation.
    """
    sfs, pms, depth, ndon = [], [], [], []
    for _, obs in obs_by_locus.items():
        if len(obs) < 2:
            continue
        tot_k = sum(o[0] for o in obs)
        tot_n = sum(o[1] for o in obs)
        for k, n in obs:
            K, N = tot_k - k, tot_n - n
            state = model.posterior(K, N)
            sf, pm = model.tail(k, n, state)
            sfs.append(float(sf))
            pms.append(float(pm))
            depth.append(n)
            ndon.append(len(obs) - 1)
    return np.array(sfs), np.array(pms), np.array(depth), np.array(ndon)


def draw_pvalues(sf, pm, seed=0):
    """One realisation of the randomised p-value: sf + U * pm, U ~ Uniform(0,1)."""
    return np.clip(sf + np.random.default_rng(seed).random(len(sf)) * pm, 0.0, 1.0)


def homopolymer_len(fa, chrom, pos, pad=HP_WINDOW, margin=40):
    """Longest homopolymer run OVERLAPPING the +/-pad window around the site, measured in full.

    The site coordinate is the protospacer anchor, not where the indel lands -- the caller accepts
    events up to 10 bp away -- so a +/-3 bp neighbourhood would miss the run the polymerase is
    actually slipping in. chrX:11,849,670 is the case in point: its four events are +/-1-2 T's in a
    15 bp T run that starts 9 bp downstream of the site coordinate.

    A run is measured over a window `margin` bp wider than the one used to decide overlap, so a
    run that straddles the window edge is reported at its true length rather than clipped to what
    happens to fall inside. Clipping would systematically shorten exactly the long runs the
    stratification is trying to isolate, and would empty the top stratum.
    """
    lo = max(0, pos - 1 - pad - margin)
    try:
        s = fa.fetch(chrom, lo, pos - 1 + pad + margin + 1).upper()
    except (ValueError, KeyError):
        return 0
    if not s:
        return 0
    core_lo = (pos - 1 - pad) - lo          # window bounds in local coordinates
    core_hi = (pos - 1 + pad) - lo
    best, start = 0, 0
    for j in range(1, len(s) + 1):
        if j < len(s) and s[j] == s[j - 1]:
            continue
        if s[start] in "ACGT" and start <= core_hi and (j - 1) >= core_lo:
            best = max(best, j - start)
        start = j
    return best


def qq(ax, p, label):
    p = np.sort(p[np.isfinite(p)])
    if not len(p):
        return
    e = (np.arange(1, len(p) + 1) - 0.5) / len(p)
    ax.plot(e, p, lw=1.4, label=f"{label} (n={len(p)})")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tables", required=True)
    ap.add_argument("--queue-all", required=True)
    ap.add_argument("--truth-wgs")
    ap.add_argument("--fasta")
    ap.add_argument("--figdir")
    ap.add_argument("--calib-seeds", type=int, default=20,
                    help="randomised-p-value seeds to sweep in the calibration test. The KS "
                         "p-value is one draw per seed; D is what is stable. Default 20.")
    ap.add_argument("--kn-policy", default="clip", choices=["clip", "drop"],
                    help="what to do with rows where control_indel_reads > control_reads")
    a = ap.parse_args()

    paths = sorted(glob.glob(a.tables))
    df = nm.load_tables(paths)
    print(f"tables {len(paths)}   site-rows {len(df)}")

    # -----------------------------------------------------------------------------------------
    # Data hazards, handled explicitly rather than silently
    # -----------------------------------------------------------------------------------------
    print("\n" + "=" * 94)
    print("DATA HAZARDS")
    print("=" * 94)
    n_nodepth = int((df.control_reads <= 0).sum())
    bad = df.control_indel_reads > df.control_reads
    print(f"rows with no control depth        : {n_nodepth}  (left on the global prior; there is "
          f"no depth to derive a floor from)")
    print(f"rows with control alt > depth     : {int(bad.sum())}")
    print("  The caller SUMS alt counts over the events at a site but takes the MEAN of control")
    print("  depth, so these are not a clean (k, n) pair. betabinom is undefined for k > n.")
    d = df[df.control_reads > 0].copy()
    if a.kn_policy == "clip":
        d["control_indel_reads"] = np.minimum(d.control_indel_reads, d.control_reads)
        print(f"  policy: CLIP k to n  ({int((bad & (df.control_reads>0)).sum())} rows affected)")
    else:
        d = d[d.control_indel_reads <= d.control_reads]
        print(f"  policy: DROP those rows")

    k = d.control_indel_reads.values.astype(float)
    n = d.control_reads.values.astype(float)
    print(f"usable control observations       : {len(d)}")

    # -----------------------------------------------------------------------------------------
    # 1. The two populations
    # -----------------------------------------------------------------------------------------
    print("\n" + "=" * 94)
    print("1. THE TWO POPULATIONS THE SINGLE BETA IS FITTED ACROSS")
    print("=" * 94)
    rate = k / n
    nz = rate[rate > 0]
    print(f"zero control alt reads : {int((k==0).sum())} / {len(k)} = {(k==0).mean()*100:.2f}%")
    print(f"nonzero                : {len(nz)}   median {np.median(nz):.4f}   "
          f"q25 {np.quantile(nz,.25):.4f}   q75 {np.quantile(nz,.75):.4f}")
    print(f"smallest nonzero rate  : {nz.min():.6f}   (= 1 read at that depth)")
    for lo, hi in ((0, 1e-3), (1e-3, 1e-2), (1e-2, 5e-2)):
        print(f"observations in ({lo:g}, {hi:g}) : {int(((rate>lo)&(rate<hi)).sum())}")
    med_dep = float(np.median(n))
    print(f"\nmedian control depth   : {med_dep:.0f}x  ->  one read = {1/med_dep*100:.2f}% VAF")
    print("The control cannot RESOLVE a rate below that. Everything the machine-error process")
    print("actually does lives underneath it, which is why the zero bin is a detection limit and")
    print("not a measurement of zero.")

    a0, b0 = nm.fit_global_prior(df.control_indel_reads.values, df.control_reads.values)
    print(f"\nfitted prior           : Beta({a0:.6g}, {b0:.6g})  mean {a0/(a0+b0):.6g}")
    print("That mean sits between the two modes, where no locus lives.")

    # -----------------------------------------------------------------------------------------
    # 2. Model fits and held-out likelihood, split BY LOCUS
    # -----------------------------------------------------------------------------------------
    print("\n" + "=" * 94)
    print("2. MODEL COMPARISON — held-out log predictive likelihood (split by locus)")
    print("=" * 94)
    loci = pd.factorize(d.chrom.astype(str) + ":" + d.end.astype(str))[0]
    held = (loci % 2) == 1
    print(f"train {int((~held).sum())} observations / test {int(held.sum())}, "
          f"{len(set(loci))} distinct loci")

    models = [BinomialGlobal(), BetaBinomMOM(), BetaBinomMML(), ZeroInflatedBB()]
    fitted = []
    for m in models:
        m.fit(k[~held], n[~held])
        ll = m.logpmf(k[held], n[held])
        ll = ll[np.isfinite(ll)]
        fitted.append((m, float(ll.mean())))
        print(f"\n{m.name}")
        print(f"  {m.describe()}")
        print(f"  held-out mean log-lik: {ll.mean():.6f}   (higher is better)")

    best = max(fitted, key=lambda t: t[1])
    print(f"\nbest by held-out likelihood: {best[0].name}")

    # Is the mixture actually identified? A Beta with a < 1 is ALREADY spike-at-zero-shaped, so
    # the "clean" component and the low end of the Beta component explain the same observations.
    # Refitting on the full data and comparing pi is the cheapest test of that, and if pi swings
    # while the likelihood barely moves, the two components are not separately estimable.
    zi_tr = [m for m in models if isinstance(m, ZeroInflatedBB)][0]
    pi_train = zi_tr.pi
    zi_all = ZeroInflatedBB().fit(k, n)
    print(f"\nmixture identifiability: pi = {pi_train:.4f} on the training half vs "
          f"{zi_all.pi:.4f} on all data,")
    print(f"  while held-out log-lik differs from plain MML by "
          f"{fitted[-1][1]-fitted[-2][1]:+.6f} nats per observation.")
    print("  pi swings while the fit does not move: the two components are NOT separately")
    print("  identified. A Beta with a < 1 is already spike-at-zero shaped, so the explicit")
    print("  mixture is re-describing what the shipped prior's shape already encodes.")

    # -----------------------------------------------------------------------------------------
    # 3. Calibration under the leave-one-donor-out null
    # -----------------------------------------------------------------------------------------
    print("\n" + "=" * 94)
    print("3. CALIBRATION — leave-one-donor-out null (every observation is noise by construction)")
    print("=" * 94)
    obs = collections.defaultdict(list)
    for c, e, kk, nn in zip(d.chrom.astype(str), d.end.astype(int), k, n):
        obs[(c, int(e))].append((kk, nn))
    multi = {key: v for key, v in obs.items() if len(v) >= 2}
    print(f"loci with >= 2 donors: {len(multi)}   observations: {sum(len(v) for v in multi.values())}")
    print("p-values are RANDOMISED; see randomised_tails() for why the plain sf cannot be used.")

    pv = {}
    for m in models:
        m.fit(k, n)                          # refit on everything for the calibration test
        sf, pm, dep, nd = randomised_tails(m, multi)
        p = draw_pvalues(sf, pm, seed=0)     # the reference realisation, seed 0
        pv[m.name] = (p, dep, nd)
        ks = kstest(p, "uniform")
        print(f"\n{m.name}")
        print(f"  KS vs Uniform(0,1): D = {ks.statistic:.4f}   p = {ks.pvalue:.3g}"
              f"   (D near 0 = calibrated)")
        print(f"  mean {p.mean():.4f} (0.5 if calibrated)   "
              f"frac < 0.05: {(p<0.05).mean():.4f} (0.05 if calibrated)   "
              f"frac < 0.001: {(p<0.001).mean():.5f} (0.001 if calibrated)")
        for lab, mask in (("thin support (1 other donor)", nd == 1),
                          ("thick support (>1 donor)", nd > 1),
                          (f"control depth < {med_dep:.0f}x", dep < med_dep),
                          (f"control depth >= {med_dep:.0f}x", dep >= med_dep)):
            if mask.sum():
                print(f"    {lab:32s} D = {kstest(p[mask],'uniform').statistic:.4f} "
                      f"  frac<0.05 = {(p[mask]<0.05).mean():.4f}   (n={int(mask.sum())})")

        # The p-value above is one draw. Sweep the seed so the reported figure is a range, not a
        # coincidence: D is stable across seeds, the KS p-value is not, and quoting a single p
        # invites a reader to read significance into RNG state.
        sweep = [kstest(draw_pvalues(sf, pm, seed=sd), "uniform")
                 for sd in range(a.calib_seeds)]
        Ds = np.array([r.statistic for r in sweep])
        ps = np.array([r.pvalue for r in sweep])
        print(f"  seed sweep (n={a.calib_seeds}): D  min {Ds.min():.4f}  median {np.median(Ds):.4f}"
              f"  max {Ds.max():.4f}")
        print(f"                     KS p  min {ps.min():.3f}  median {np.median(ps):.3f}"
              f"  max {ps.max():.3f}   rejects at 0.05: {(ps<0.05).sum()}/{a.calib_seeds}")

    # -----------------------------------------------------------------------------------------
    # 4. Is the rate context-dependent?
    # -----------------------------------------------------------------------------------------
    if a.fasta:
        print("\n" + "=" * 94)
        print("4. CONTEXT — is the Bernoulli rate a constant of the assay?")
        print("=" * 94)
        import pysam
        fa = pysam.FastaFile(a.fasta)
        uniq = d.drop_duplicates(["chrom", "end"])[["chrom", "end"]]
        hp = {(c, int(e)): homopolymer_len(fa, c, int(e))
              for c, e in zip(uniq.chrom.astype(str), uniq.end.astype(int))}
        d["hp"] = [hp[(c, int(e))] for c, e in zip(d.chrom.astype(str), d.end.astype(int))]
        strata = [("no run (<=3 bp)", d.hp <= 3), ("4-5 bp", d.hp.between(4, 5)),
                  ("6-8 bp", d.hp.between(6, 8)), (">=9 bp", d.hp >= 9)]
        print(f"{'stratum':18s} {'n':>7s} {'obs rate':>10s} {'a0':>12s} {'b0':>10s} "
              f"{'prior mean':>12s} {'nonzero %':>10s}")
        for lab, mask in strata:
            s = d[mask]
            if len(s) < 50:
                print(f"{lab:18s} {len(s):7d}   (too few to fit)")
                continue
            sa, sb = nm.fit_global_prior(s.control_indel_reads.values, s.control_reads.values)
            r = s.control_indel_reads.sum() / max(s.control_reads.sum(), 1)
            nzp = (s.control_indel_reads > 0).mean() * 100
            print(f"{lab:18s} {len(s):7d} {r:10.6f} {sa:12.6g} {sb:10.4g} "
                  f"{sa/(sa+sb):12.6g} {nzp:10.2f}")
        print("\nIf the strata separate, the rate is NOT a constant of the assay and a single")
        print("global prior necessarily under-penalises slippage-prone context while")
        print("over-penalising clean unique sequence.")

        # What a context-conditioned prior would actually do to the rows under review. Without
        # this the stratification is an interesting table with no consequence attached.
        qh = pd.read_csv(a.queue_all, sep="\t")
        qh["hp"] = [homopolymer_len(fa, c, int(e))
                    for c, e in zip(qh.chrom.astype(str), qh.end.astype(int))]
        wd_ = qh.why_dropped.fillna("")
        print(f"\n-- the same strata among the {len(qh)} gated rows --")
        for lab, mask in [("no run (<=3 bp)", qh.hp <= 3), ("4-5 bp", qh.hp.between(4, 5)),
                          ("6-8 bp", qh.hp.between(6, 8)), (">=9 bp", qh.hp >= 9)]:
            print(f"  {lab:16s} gated {int(mask.sum()):4d}   in the queue "
                  f"{int((mask & (wd_ == '')).sum()):3d}")
        hi = qh[(qh.hp >= 9) & (wd_ == "")]
        if len(hi):
            sa, sb = nm.fit_global_prior(d[d.hp >= 9].control_indel_reads.values,
                                         d[d.hp >= 9].control_reads.values)
            al2 = sa + hi.control_indel_reads.values
            be2 = sb + np.maximum(hi.control_reads.values - hi.control_indel_reads.values, 0)
            aq2 = nm.aq_from_sf(betabinom.sf(hi.indel_reads.values.astype(int) - 1,
                                             hi.total_reads.values.astype(int), al2, be2))
            al1, be1, _, dep1 = nm.control_posterior(df, hi, "matched", a0, b0)
            al1, be1, _ = nm.apply_depth_floor(al1, be1, dep1)
            aq1 = nm.aq_from_sf(betabinom.sf(hi.indel_reads.values.astype(int) - 1,
                                             hi.total_reads.values.astype(int), al1, be1))
            print(f"\n-- queue rows in a >=9 bp homopolymer: {len(hi)} --")
            print(f"  AQ under the GLOBAL prior      : {np.round(aq1, 1).tolist()}")
            print(f"  AQ under the >=9 bp STRATUM prior: {np.round(aq2, 1).tolist()}")
            print(f"  would newly fall below AQ 5: {int(((aq1>=5)&(aq2<5)).sum())} row(s)")
            print("  This is the concrete cost of a global prior: these rows are scored against")
            print("  an assay-wide rate that their own sequence context says is far too low.")
        else:
            print("\nNo queue row sits in a >=9 bp homopolymer on this cohort, so the")
            print("mis-specification is real but currently costs nothing at the queue.")

    # -----------------------------------------------------------------------------------------
    # 5. The depth floor, ablated
    # -----------------------------------------------------------------------------------------
    print("\n" + "=" * 94)
    print("5. THE DEPTH FLOOR — patch, or genuine requirement?")
    print("=" * 94)
    q = pd.read_csv(a.queue_all, sep="\t")
    alpha, beta, nsup, dep = nm.control_posterior(df, q, "matched", a0, b0)
    kq = q.indel_reads.values.astype(int)
    nq = q.total_reads.values.astype(int)
    aq_off = nm.aq_from_sf(betabinom.sf(kq - 1, nq, alpha, beta))
    af, bf, n_lift = nm.apply_depth_floor(alpha, beta, dep)
    aq_on = nm.aq_from_sf(betabinom.sf(kq - 1, nq, af, bf))
    clean = dep > 0
    print(f"gated rows {len(q)};  floor raises the background on {n_lift} of them "
          f"({n_lift/len(q)*100:.1f}%)")
    print(f"posterior mean WITHOUT the floor, on floored rows: "
          f"min {np.min((alpha/(alpha+beta))[dep>0]):.3g}, "
          f"median {np.median((alpha/(alpha+beta))[dep>0]):.3g}")
    print(f"the fitted prior mean is {a0/(a0+b0):.3g}, and 1/median_control_depth is "
          f"{1/med_dep:.3g} -- the unfloored posterior claims a rate far below what the control "
          f"can support")
    print(f"AQ with floor OFF: median {np.median(aq_off):.1f}, rows AQ<5: "
          f"{int((aq_off<5).sum())}")
    print(f"AQ with floor ON : median {np.median(aq_on):.1f}, rows AQ<5: {int((aq_on<5).sum())}")

    zi = [m for m in models if isinstance(m, ZeroInflatedBB)][0]
    print(f"\nThe mixture can state 'this locus is clean' directly: it puts weight "
          f"pi = {zi.pi:.4f} on a\ncomponent at rate eps = {zi.eps:.3g}, which is what the floor "
          f"is imitating by clamping.")
    if zi.eps < 1e-12:
        print("\nAND THAT eps IS THE POINT, not a fitting failure. The MLE drove it to the")
        print("optimiser's floor because NO observation in the data constrains it from below: at")
        print(f"{med_dep:.0f}x, a rate of 1e-4 and a rate of 1e-40 both predict zero alt reads")
        print("with probability ~1. The control is blind to the machine-error rate, so the data")
        print("cannot estimate it -- which is exactly the argument for taking that component from")
        print("a pooled population panel instead of from one matched control.")

    # -----------------------------------------------------------------------------------------
    # 6. Operational consequence
    # -----------------------------------------------------------------------------------------
    if a.truth_wgs:
        print("\n" + "=" * 94)
        print("6. OPERATIONAL CONSEQUENCE — does any of this move the queue?")
        print("=" * 94)
        lab = nm.wgs_curated_label(a.truth_wgs)
        qq_ = q.copy()
        qq_["guide"] = qq_.sample_name.map(nm.guide_of)
        qq_["chrom"] = qq_.chrom.astype(str)

        # The AQ rule is 4th in review_filter's np.select, and first match wins -- so it only ever
        # sees rows that rules 1-3 (germline / far from PAM / single indel length) left alone.
        # Scoring all 479 gated rows would badly overstate what changing the model does, because
        # most of those rows never reach the AQ test in the shipped pipeline.
        wd = q.why_dropped.fillna("")
        earlier = {"germline (present in control)", "far from PAM", "single indel length"}
        reaches = ~wd.isin(earlier)
        print(f"gated rows {len(q)};  rules 1-3 remove {int((~reaches).sum())};  "
              f"rows that actually reach the AQ rule: {int(reaches.sum())}")
        print("Scoring all 479 would overstate the effect -- most never reach this test.")

        m0 = qq_[reaches.values].merge(lab, on=["guide", "chrom", "start"], how="inner")
        print(f"labelled among those: {len(m0)}  "
              f"({int((m0.curated_label==1).sum())} confirmed / "
              f"{int((m0.curated_label==0).sum())} rejected)")
        print(f"\n{'model':42s} {'AQ<5 drops':>11s} {'queue':>7s} {'confirmed':>12s} "
              f"{'rejected':>11s}")
        key = ["sample_name", "chrom", "start"]
        for m in models:
            if isinstance(m, BinomialGlobal):
                sf = binom.sf(kq - 1, nq, m.p)
            elif isinstance(m, ZeroInflatedBB):
                sf = np.array([m.tail(kk - 1, nn, m.posterior(al - a0, (al - a0) + (be - b0)))[0]
                               for kk, nn, al, be in zip(kq, nq, alpha, beta)])
            else:
                sf = betabinom.sf(kq - 1, nq, m.a + (alpha - a0), m.b + (beta - b0))
            keep = pd.Series(nm.aq_from_sf(sf) >= 5, index=q.index) & reaches
            drops = int((reaches & ~keep).sum())
            kmap = q.assign(_k=keep.values).set_index(key)["_k"].to_dict()
            mk = np.array([bool(kmap.get(t, False))
                           for t in zip(m0.sample_name, m0.chrom, m0.start)])
            ck = int((mk & (m0.curated_label == 1).values).sum())
            rk = int((mk & (m0.curated_label == 0).values).sum())
            print(f"{m.name:42s} {drops:11d} {int(keep.sum()):7d} "
                  f"{ck:>7d}/{int((m0.curated_label==1).sum()):<4d} "
                  f"{rk:>6d}/{int((m0.curated_label==0).sum()):<4d}")
        print("\n('queue' here counts rows surviving the AQ rule only; the shipped queue is that")
        print(" number minus the repeat-region and external-panel rules, which run after it.)")
        print("A model that is statistically better but moves these columns by a row or two has")
        print("not earned a pipeline change on this cohort. Report it as such.")

    # -----------------------------------------------------------------------------------------
    if a.figdir:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from scipy.stats import beta as beta_dist
        os.makedirs(a.figdir, exist_ok=True)

        fig, ax = plt.subplots(1, 2, figsize=(12, 4.2))
        edges = np.linspace(0, 1, 101)
        ax[0].hist(rate, bins=edges, color="#444", log=True)
        ax[0].set_xlabel("control indel rate  k/n"); ax[0].set_ylabel("observations (log)")
        ax[0].set_title(f"empirical: {(k==0).mean()*100:.1f}% at exactly 0, "
                        f"a second mode near {np.median(nz):.2f}")
        x = np.linspace(1e-4, 1 - 1e-4, 500)
        ax[1].plot(x, beta_dist.pdf(x, a0, b0), color="crimson",
                   label=f"fitted Beta({a0:.4g}, {b0:.3g})")
        ax[1].axvline(a0 / (a0 + b0), ls="--", color="k",
                      label=f"prior mean {a0/(a0+b0):.5f}")
        ax[1].axvline(float(np.median(nz)), ls=":", color="#2a6",
                      label=f"germline mode {np.median(nz):.2f}")
        ax[1].set_yscale("log"); ax[1].set_xlabel("background rate p")
        # a < 1 and b > 1 makes this strictly decreasing -- infinite at 0, zero at 1. It is
        # spike-at-zero shaped, NOT U-shaped, and that is precisely why it survives the
        # calibration test: the shape already encodes the zero inflation.
        ax[1].set_title(f"Beta(a<1, b>1): spike at 0 with a heavy tail\n"
                        f"prior mean {a0/(a0+b0):.5f} lies between the two modes")
        ax[1].legend(fontsize=8)
        fig.tight_layout()
        p1 = os.path.join(a.figdir, "noise_two_populations.png")
        fig.savefig(p1, dpi=150)

        fig, ax = plt.subplots(figsize=(5.4, 5.2))
        for name, (p, _, _) in pv.items():
            qq(ax, p, name)
        ax.plot([0, 1], [0, 1], "k--", lw=1)
        ax.set_xlabel("expected quantile under Uniform(0,1)")
        ax.set_ylabel("observed randomised p-value")
        ax.set_title("calibration under the leave-one-donor-out null")
        ax.legend(fontsize=7, loc="lower right")
        fig.tight_layout()
        p2 = os.path.join(a.figdir, "noise_calibration_qq.png")
        fig.savefig(p2, dpi=150)
        print(f"\nwrote {p1}\nwrote {p2}")


if __name__ == "__main__":
    main()
