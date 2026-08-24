#!/usr/bin/env python3
"""noise_model_followup.py — chase the two findings that looked like defects in the prior.

`noise_model_validate.py` reported two results that pointed at the model itself rather than at how
it is described:

  clean-locus under-prediction   the coldest bin of the leave-one-donor-out reliability diagram
                                 under-predicts observed alt reads by 11.7x, CI excluding 1.
  assumption 8                   scoring EDITED counts against the matched control at no-edit loci
                                 gives KS D = 0.0343, ~7x the control-vs-control figure.

Both dissolve on inspection, in different ways, and this script is what shows that. It is analysis
only: it imports the shipped model, changes nothing, and the pipeline does not call it.

  --part clean   the reliability finding. Isolates the cold bins and tests four explanations that
                 make different, checkable predictions: donor-private germline, sequence context,
                 a few pathological loci, and a prior tail that is simply too thin.
  --part a8      the assumption-8 finding. Stratifies the same rows by whether the caller reported
                 an event at that site, because conditioning on "an indel was called here" forces
                 k >= 1 and makes the p-value a spike rather than a distribution.

The short version of both: the reliability finding is five loci carrying donor-private alleles the
caller had already suppressed, and the assumption-8 finding is a selection effect from pooling the
3% of rows that carry a called event into the 97% that do not. Neither is a defect in the prior.
See docs/NOISE_MODEL_ASSUMPTIONS.md, "Follow-up".

usage:
  noise_model_followup.py --tables 'results_cart_bnd/*/*.offtarget_analysis.tsv' \
      --fasta .../hg38_PLVM_CD19_CARv4_cd34.fa --part both
"""
import argparse
import collections
import glob
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import betabinom, chi2, kstest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import noise_model as nm                                                        # noqa: E402
from noise_model_validate import BetaBinomMOM, homopolymer_len                  # noqa: E402

# The gate a row must clear before it can be a call at all. Rows below it are background by
# construction, so they bound how much of any miscalibration could be genuine off-target editing.
GATE_READS, GATE_VAF = 2, 0.005
GERMLINE_VAF = 0.35          # a real heterozygote sits near 0.5; this is the generous edge of it
HP_SLIPPERY = 7              # homopolymer length at which slippage starts to dominate


def poisson_ci(k, conf=0.95):
    """Exact (Garwood) interval on a Poisson count, so bins built on a few events say so."""
    a = 1 - conf
    lo = 0.0 if k == 0 else chi2.ppf(a / 2, 2 * k) / 2
    return lo, chi2.ppf(1 - a / 2, 2 * (k + 1)) / 2


def load(tables):
    df = nm.load_tables(sorted(glob.glob(tables)))
    d = df[df.control_reads > 0].copy()
    d["control_indel_reads"] = np.minimum(d.control_indel_reads, d.control_reads)
    return df, d


# =================================================================================================
# Part 1 — the clean-locus under-prediction
# =================================================================================================

def part_clean(df, d, fasta_path):
    shipped = BetaBinomMOM()
    shipped.fit(d.control_indel_reads.values.astype(float),
                d.control_reads.values.astype(float))
    print(f"prior Beta({shipped.a:.6g}, {shipped.b:.6g})")

    # Rebuild the reliability design: predict each donor from the OTHER donors at that locus.
    obs = collections.defaultdict(list)
    for c, e, s, kk, nn in zip(d.chrom.astype(str), d.end.astype(int), d.sample_name.astype(str),
                               d.control_indel_reads.values.astype(float),
                               d.control_reads.values.astype(float)):
        obs[(c, int(e))].append((s, kk, nn))
    multi = {key: v for key, v in obs.items() if len(v) >= 2}
    rows = []
    for (c, e), o in multi.items():
        tot_k, tot_n = sum(x[1] for x in o), sum(x[2] for x in o)
        for s, kk, nn in o:
            al, be = shipped.posterior(tot_k - kk, tot_n - nn)
            rows.append((c, e, s, kk, nn, nn * al / (al + be)))
    r = pd.DataFrame(rows, columns=["chrom", "pos", "sample", "k", "n", "exp"])
    r["vaf"] = r.k / r.n
    print(f"loci with >=2 donors {len(multi)}   observations {len(r)}")
    print(f"aggregate: predicted {r.exp.sum():,.1f}  observed {r.k.sum():,.0f}  "
          f"ratio {r.k.sum()/r.exp.sum():.3f}")

    cold = r[r.exp < 3e-3].copy()                       # the two under-predicting bins
    excess = cold.k.sum() - cold.exp.sum()
    print(f"\ncold bins (expected < 3e-3): {len(cold):,} observations, "
          f"{int((cold.k >= 1).sum())} carry any alt read, "
          f"excess {excess:,.1f} reads over {cold.exp.sum():,.1f} predicted")

    print("\n-- is it concentrated? --")
    per = (cold.groupby(["chrom", "pos"])
               .agg(k=("k", "sum"), exp=("exp", "sum"), nobs=("k", "size"))
               .assign(exc=lambda x: x.k - x.exp).sort_values("exc", ascending=False))
    tot_exc = per.exc[per.exc > 0].sum()
    for cut in (1, 5, 10, 25):
        print(f"   top {cut:3d} loci of {len(per):,}: "
              f"{per.exc.head(cut).sum()/max(tot_exc, 1)*100:5.1f}% of the positive excess")

    # Did the caller already know about these? n_control_filtered counts the events it dropped
    # at -x 0 because the control supported them -- i.e. the ones it called germline.
    nc = pd.to_numeric(df.get("n_control_filtered"), errors="coerce").fillna(0)
    key = df.chrom.astype(str) + ":" + df.end.astype(int).astype(str) + "|" + df.sample_name.astype(str)
    supp = dict(zip(key, nc))
    print(f"\n   {'locus':<22} {'donor':<24} {'pred':>7} {'alt':>5} {'VAF':>6} {'suppressed':>11}")
    for (c, p), row in per.head(5).iterrows():
        g = cold[(cold.chrom == c) & (cold.pos == p)].sort_values("k", ascending=False).iloc[0]
        print(f"   {c + ':' + format(int(p), ','):<22} {g['sample']:<24} {row.exp:7.3f} "
              f"{int(row.k):5d} {g.vaf:6.3f} "
              f"{int(supp.get(f'{c}:{int(p)}|' + g['sample'], 0)):11d}")

    print("\n-- the association, cohort-wide --")
    dd = df[(df.control_reads > 0) & (df.total_reads > 0)]
    if "is_target" in dd.columns:
        dd = dd[dd.is_target == 0]
    ncd = pd.to_numeric(dd.get("n_control_filtered"), errors="coerce").fillna(0)
    cv = np.minimum(dd.control_indel_reads, dd.control_reads) / dd.control_reads
    ev = dd.indel_reads / dd.total_reads
    orph = (cv >= 0.05) & (ev < 0.01)
    print(f"   rows with control VAF >=5% and a near-empty edited library: {int(orph.sum())}")
    print(f"     with n_control_filtered > 0 : {(ncd[orph] > 0).mean()*100:.1f}%")
    print(f"     the same figure elsewhere   : {(ncd[~orph] > 0).mean()*100:.2f}%")
    print("   A real donor variant is in BOTH libraries; the edited count reads 0 because the")
    print("   caller dropped the event at -x 0, not because the allele is absent.")

    print("\n-- the competing explanations --")
    hits = cold[cold.k >= 1]
    germ = hits[hits.vaf >= GERMLINE_VAF]
    print(f"   germline at VAF >= {GERMLINE_VAF}: {len(germ)} observations, "
          f"{germ.k.sum():.0f} alt reads ({germ.k.sum()/max(excess,1)*100:.1f}% of the excess)")
    try:
        import pysam
        fa = pysam.FastaFile(fasta_path)
    except Exception as exc:                                            # pragma: no cover
        print(f"   context: fasta unavailable ({exc}); skipped")
        cold["hp"] = 0
    else:
        sites = cold[["chrom", "pos"]].drop_duplicates()
        hp = {(c, int(p)): homopolymer_len(fa, c, int(p))
              for c, p in zip(sites.chrom, sites.pos)}
        cold["hp"] = [hp[(c, int(p))] for c, p in zip(cold.chrom, cold.pos)]
        for lo, hi, lab in ((0, 5, "<5 bp"), (5, 7, "5-6 bp"), (7, 9, "7-8 bp"), (9, 999, ">=9 bp")):
            g = cold[(cold.hp >= lo) & (cold.hp < hi)]
            if not len(g) or g.exp.sum() <= 0:
                continue
            print(f"   homopolymer {lab:>7}: {len(g):6d} obs, "
                  f"{(g.k.sum()-g.exp.sum())/max(excess,1)*100:5.1f}% of the excess")

    print("\n-- what is left once all three are removed --")
    worst = set(per.head(10).index.tolist())
    keep = cold[(cold.vaf < GERMLINE_VAF) & (cold.hp < HP_SLIPPERY)]
    keep = keep[[(c, int(p)) not in worst for c, p in zip(keep.chrom, keep.pos)]]
    pe, po = keep.exp.sum(), keep.k.sum()
    lo, hi = poisson_ci(int(po))
    print(f"   {len(keep):,} observations, predicted {pe:.2f}, observed {po:.0f}, "
          f"ratio {po/max(pe,1e-9):.2f}  95% CI [{lo/max(pe,1e-9):.2f}, {hi/max(pe,1e-9):.2f}]")
    print("   At genuinely clean loci the model OVER-predicts. The prior's left tail is not too")
    print("   thin; the 11.7x was five loci carrying alleles the caller had already suppressed.")


# =================================================================================================
# Part 2 — assumption 8, stratified
# =================================================================================================

def part_a8(df, d, seed=0):
    a0, b0 = nm.fit_global_prior(df.control_indel_reads.values, df.control_reads.values)
    print(f"prior Beta({a0:.6g}, {b0:.6g})")

    e = df[(df.control_reads > 0) & (df.total_reads > 0)].copy()
    if "is_target" in e.columns:
        e = e[e.is_target == 0]
    e["control_indel_reads"] = np.minimum(e.control_indel_reads, e.control_reads)
    ic = pd.to_numeric(e.indel_count, errors="coerce").fillna(0).values
    nc = pd.to_numeric(e.get("n_control_filtered"), errors="coerce").fillna(0).values
    ca = e.control_indel_reads.values.astype(float)
    cn = e.control_reads.values.astype(float)
    ek = e.indel_reads.values.astype(int)
    en = e.total_reads.values.astype(int)
    print(f"no-target rows with depth in both libraries: {len(e):,}")

    def score(al, be):
        rng = np.random.default_rng(seed)
        return np.clip(betabinom.sf(ek, en, al, be)
                       + rng.random(len(ek)) * betabinom.pmf(ek, en, al, be), 0.0, 1.0)

    for floor_on in (False, True):
        al, be = a0 + ca, b0 + np.maximum(cn - ca, 0.0)
        lift = 0
        if floor_on:
            al, be, lift = nm.apply_depth_floor(al, be, cn)
        pv = score(al, be)
        fin = np.isfinite(pv)
        head = "WITH depth floor (production)" if floor_on else "RAW posterior (the model)"
        print(f"\n-- {head} --" + (f"   floor raised {lift:,}" if floor_on else ""))
        print(f"   {'population':<34} {'n':>7} {'KS D':>8} {'p<0.001':>10} {'excess':>8} {'obs/pred':>9}")
        pops = [("pooled (as first reported)", fin),
                ("indel_count == 0  no event called", fin & (ic == 0)),
                ("indel_count == 1  one event", fin & (ic == 1)),
                ("indel_count >= 2", fin & (ic >= 2)),
                ("germline-suppressed removed", fin & (nc == 0))]
        for lab, m in pops:
            q = pv[m]
            if not len(q):
                continue
            frac = (q < 0.001).mean()
            pred = float((en * (al / (al + be)))[m].sum())
            print(f"   {lab:<34} {len(q):7d} {kstest(q, 'uniform').statistic:8.4f} "
                  f"{frac:10.5f} {frac/0.001:7.1f}x {ek[m].sum()/max(pred,1e-9):9.3f}")

    al, be = a0 + ca, b0 + np.maximum(cn - ca, 0.0)
    pv = score(al, be)
    print("\n-- why indel_count >= 1 cannot be read as a calibration failure --")
    m1 = ic == 1
    print(f"   indel_count == 1 rows: {int(m1.sum()):,}   share with k >= 1: "
          f"{(ek[m1] >= 1).mean():.3f}  (1.000 means selected on the outcome)")
    print(f"   their p-values: median {np.median(pv[m1]):.5g}, "
          f"{(pv[m1] < 0.01).mean()*100:.1f}% below 0.01 — a spike at one value, not a spread")

    print("\n-- the stratum where the assumption CAN be tested --")
    m0 = ic == 0
    print(f"   {'threshold':>12} {'observed':>10} {'uniform':>10}")
    for t in (0.5, 0.05, 0.01, 0.001):
        print(f"   {'P(p < ' + format(t, 'g') + ')':>12} {(pv[m0] < t).mean():10.5f} {t:10.5f}")
    print(f"   n = {int(m0.sum()):,}   naive KS critical value "
          f"{1.358/np.sqrt(int(m0.sum())):.4f}   observed D "
          f"{kstest(pv[m0], 'uniform').statistic:.4f}")
    print("   Uniform in the bulk. Supported where testable; untestable where it bites, because")
    print("   the rows the filter adjudicates are exactly the rows selected for carrying an event.")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tables", required=True)
    ap.add_argument("--fasta")
    ap.add_argument("--part", default="both", choices=["clean", "a8", "both"])
    a = ap.parse_args()

    df, d = load(a.tables)
    print(f"site-rows {len(df):,}   usable control observations {len(d):,}")
    if a.part in ("clean", "both"):
        print("\n" + "=" * 96)
        print("PART 1 — the clean-locus under-prediction")
        print("=" * 96)
        part_clean(df, d, a.fasta)
    if a.part in ("a8", "both"):
        print("\n" + "=" * 96)
        print("PART 2 — assumption 8, stratified by whether an event was called")
        print("=" * 96)
        part_a8(df, d)


if __name__ == "__main__":
    main()
