#!/usr/bin/env python3
"""Build a panel of normals (PoN) of bad off-target hotspot sites.

A predicted off-target site that carries indel reads in a sample that was never edited cannot
be a CRISPR edit. Such sites are germline variants or recurrent alignment artifacts -- the
off-target panel is built from homology, so it is deliberately enriched for repetitive and
paralogous sequence, exactly where aligners invent indels.

Input is `find_edited_reads.py` output (*.offtarget_analysis.tsv) produced with an UNEDITED
sample as the query. Output is a small TSV of blacklisted sites consumed by review_filter.py.

This is the single-guide-safe form of "this site lights up regardless of which guide was used".
The cross-guide recurrence rule needs a multi-guide cohort to work at all; a PoN does not, and
on the CAR-T cohort it is strictly better (queue 63 vs 65, 2 false positives vs 4).

usage:
  build_offtarget_pon.py unedited1.tsv unedited2.tsv ... -o assets/offtarget_pon.tsv
"""
import argparse
import os
import sys

import pandas as pd

MIN_READS = 3
MIN_DONORS = 1
# A read-count threshold alone is not depth-robust. At WGS depth (~60x) 3 reads is ~5% VAF, a
# real germline/artifact signal. At ECS depth (~2000x) 3 reads is ~0.15% -- sequencing noise --
# and blacklists half the panel. The VAF floor keeps the same rule meaningful in both arms.
MIN_VAF = 0.02


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("inputs", nargs="+",
                    help="*.offtarget_analysis.tsv scored with an UNEDITED sample as query")
    ap.add_argument("-o", "--out", required=True, help="PoN TSV to write")
    ap.add_argument("--min-reads", type=int, default=MIN_READS,
                    help="indel reads in an unedited sample for the site to count (default 3)")
    ap.add_argument("--min-donors", type=int, default=MIN_DONORS,
                    help="number of distinct unedited samples that must show it (default 1)")
    ap.add_argument("--min-vaf", type=float, default=MIN_VAF,
                    help="indel fraction in an unedited sample for the site to count "
                         "(default 0.02; keeps the rule meaningful at ECS depth)")
    ap.add_argument("--merge", help="an existing PoN to carry forward; sites are keyed by "
                                    "position, so this accumulates recurrent bad regions across "
                                    "runs without discarding either panel's coverage")
    a = ap.parse_args()

    frames = []
    for p in a.inputs:
        d = pd.read_csv(p, sep="\t")
        d["donor"] = os.path.basename(p).split(".offtarget")[0]
        frames.append(d)
    n = pd.concat(frames, ignore_index=True)

    donors = n.donor.nunique()
    if donors < 2:
        print(f"WARNING: only {donors} unedited sample(s) supplied. A PoN from a single donor "
              f"cannot separate that donor's private germline from recurrent artifact; the "
              f"result will still work but will be donor-biased.", file=sys.stderr)

    # Emit EVERY scored site, not just the bad ones, with a `blacklisted` flag.
    #
    # The blacklist alone is unusable as a safety check: a PoN built for one guide panel shares
    # almost no coordinates with a different guide's panel, so applying the wrong PoN looks
    # exactly like applying a clean one -- zero hits, no error. Carrying the full scored universe
    # lets review_filter.py measure what fraction of the sites it is filtering this PoN actually
    # covers, and refuse to pretend rule 4 ran when it did not.
    # A site counts against a donor only if that donor shows both enough reads and enough VAF.
    n["is_bad"] = (n.indel_reads >= a.min_reads) & (n.indel_fraction >= a.min_vaf)
    universe = (n.assign(bad_donor=n.donor.where(n.is_bad))
                  .groupby(["chrom", "end"])
                  .agg(n_donors=("bad_donor", "nunique"),
                       max_indel_reads=("indel_reads", "max"),
                       max_indel_fraction=("indel_fraction", "max"))
                  .reset_index()
                  .rename(columns={"end": "pos"}))
    universe["blacklisted"] = (universe.n_donors >= a.min_donors).astype(int)
    n_new = len(universe)
    if a.merge and os.path.exists(a.merge):
        prev = pd.read_csv(a.merge, sep="\t")
        if "blacklisted" not in prev.columns:       # tolerate a pre-universe PoN
            prev["blacklisted"] = 1
        universe = pd.concat([prev, universe], ignore_index=True)
        # A site is blacklisted if EITHER panel found it bad; keep the strongest evidence seen.
        universe = (universe.groupby(["chrom", "pos"])
                            .agg(n_donors=("n_donors", "max"),
                                 max_indel_reads=("max_indel_reads", "max"),
                                 max_indel_fraction=("max_indel_fraction", "max"),
                                 blacklisted=("blacklisted", "max"))
                            .reset_index())
        print(f"merged with {a.merge}: {len(prev)} prior sites + {n_new} new "
              f"-> {len(universe)} total")

    universe = universe.sort_values(["chrom", "pos"])
    universe.to_csv(a.out, sep="\t", index=False)

    n_bad = int(universe.blacklisted.sum())
    print(f"unedited samples  : {donors} ({', '.join(sorted(n.donor.unique()))})")
    print(f"site-rows scored  : {len(n)}")
    print(f"panel sites       : {len(universe)}   (the coverage universe)")
    print(f"blacklisted sites : {n_bad} "
          f"(>={a.min_reads} reads AND >={a.min_vaf:.0%} VAF in >={a.min_donors} unedited sample)  -> {a.out}")
    print()
    print("This PoN is specific to the guide panel it was built from. To use the filter with a "
          "new\nguide, score your unedited sample(s) against THAT guide's target file and "
          "rebuild.")


if __name__ == "__main__":
    main()
