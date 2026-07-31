#!/usr/bin/env python3
"""
indel_length_diversity.py — measure the "varied indel lengths at one cut" signature.

The heuristic this quantifies
-----------------------------
From manual IGV review: **if every supporting read shows the exact same indel length, it is
probably a sequencing/alignment artifact.** A genuine CRISPR pileup carries indels of *different*
sizes at *slightly* different places, all overlapping the one position Cas9 cut. NHEJ repair is
stochastic, so real editing is heterogeneous; a systematic artifact is monotonous.

The shape model has no such feature — it has `modal_len` (the mode) and nothing about the spread.
Note this is NOT the same thing as `conc_ratio`/`pos_conc`, which reward *positional* concordance.
Both are wanted, and they pull in opposite directions:

    position should be CONCORDANT   (one cut site)
    length   should be DIVERSE      (stochastic repair)

Measured on the AAVS1 cohort (from the WGS tagged BAMs) — and the honest answer is that on
THIS cohort it does not yet separate:

    raw distinct-length count   TP (n=2): 19, 16      TN (n=26): median 2, max 14
    depth-controlled at 30      TP:        8,  9      TN:        9, 5, 2

The raw count looks decisive, but it is confounded with depth: the TP sites carry 119-149
edit-supporting reads while most TN sites carry 1-17. Only three TN sites have >=30 edit reads,
and once every site is rarefied to the same 30 reads the best TN (chr6:36797655) scores 9 --
matching the two real edits at 8 and 9. **No separation survives the depth control.** This is the
same trap that removed `n_distinct_pos` from the old feature set (`features.py:12`, "it grows").

Which does not mean the heuristic is wrong. The two TN sites with the highest depth-controlled
diversity, chr6:36797655 and chr7:26760750, are also the two with the highest ECS VAF (0.248 and
0.302). They may be genuine off-target edits that curation never confirmed rather than artifacts,
in which case the feature is behaving correctly and the labels are incomplete. chr6:36797655 is
the site to open first in IGV: high length diversity and strong ECS support say "real", while
sitting 41 bp from the predicted cut says "artifact".

**Do not use normalised Shannon entropy for this.** It divides by log(n_distinct), which cancels
the very signal we want: measured TP normalised entropy is 0.55-0.59 while several TN sites reach
1.000 — the ranking inverts. Report the raw distinct-length count together with the edit-read
depth and let the tree model learn the interaction; gradient boosting handles that natively and
needs no hand-chosen normalisation.

Usage
-----
  indel_length_diversity.py --outdir results_offtarget_aavs1_igv/offtarget [--out div.csv]
"""
import argparse
import collections
import math
import os
import random
import re
import statistics
import sys

import pandas as pd
import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from validate_recall_aavs1 import load_truth, DEFAULT_TP, DEFAULT_TN, SLACK, _near

TAG_LEN = re.compile(r"Edited_(Deletion|Insertion)_(\d+)bp")
WINDOW = 60


def edit_lengths(bam, chrom, pos, window=WINDOW, tag="XC"):
    """Indel lengths of the edit-supporting reads around a site, from their XC tags."""
    out = []
    try:
        it = bam.fetch(chrom, max(0, pos - window), pos + window)
    except (ValueError, KeyError):
        return out
    for r in it:
        try:
            m = TAG_LEN.match(r.get_tag(tag) or "")
        except KeyError:
            continue
        if m:
            out.append(int(m.group(2)))
    return out


def summarise(lengths):
    if not lengths:
        return None
    c = collections.Counter(lengths)
    n = sum(c.values())
    return {"n_edit_reads": n,
            "n_distinct_len": len(c),
            "modal_len_frac": round(max(c.values()) / n, 3),
            # unnormalised Shannon entropy: unlike the normalised form it keeps the
            # "how many different lengths" signal instead of dividing it out
            "len_entropy": round(-sum((v / n) * math.log(v / n) for v in c.values()), 3)}


def rarefy(lengths, k, draws=200, seed=0):
    """Distinct lengths at a FIXED subsample size — the depth-controlled comparison."""
    if len(lengths) < k:
        return None
    rng = random.Random(seed)
    return statistics.median(len(set(rng.sample(lengths, k))) for _ in range(draws))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", required=True, help="the run's offtarget/ directory")
    ap.add_argument("--tp", default=DEFAULT_TP)
    ap.add_argument("--tn", default=DEFAULT_TN)
    ap.add_argument("--rarefy-to", type=int, default=30,
                    help="depth-controlled distinct-length count at this many edit reads")
    ap.add_argument("--out", help="write the per-site table here")
    args = ap.parse_args()

    tp_set, tn_set = load_truth(args.tp, args.tn)
    t = pd.read_csv(os.path.join(args.outdir, "training.tsv"), sep="\t", low_memory=False)
    t["called"] = t["verdict"].astype(str).str.contains("LIKELY EDIT", na=False)
    t["curated"] = ["TP" if _near(c, p, tp_set, SLACK) else
                    ("TN" if _near(c, p, tn_set, SLACK) else "uncurated")
                    for c, p in zip(t["chrom"], t["start"])]

    tag_dir = os.path.join(args.outdir, "wgs_tagged")
    if not os.path.isdir(tag_dir):
        sys.exit(f"ERROR: {tag_dir} not found — rerun with --offtarget_wgs_tagged_bam true")
    cache = {}
    rows = []
    for _, r in t[t["called"] | (t["curated"] == "TP")].iterrows():
        s = r["sample"]
        if s not in cache:
            p = os.path.join(tag_dir, f"{s}.wgs_tagged.bam")
            cache[s] = pysam.AlignmentFile(p) if os.path.exists(p) else None
        if cache[s] is None:
            continue
        pos = int(r["modal_pos"]) if pd.notna(r.get("modal_pos")) else int(r["start"])
        lengths = edit_lengths(cache[s], str(r["chrom"]), pos)
        summ = summarise(lengths)
        if summ is None:
            continue
        rows.append({"curated": r["curated"], "site": f"{r['chrom']}:{int(r['start'])}",
                     "cut_dist": r.get("cut_dist"), "ecs_if": r.get("ecs_if"),
                     **summ,
                     f"distinct_at_{args.rarefy_to}": rarefy(lengths, args.rarefy_to)})

    df = pd.DataFrame(rows).sort_values(["curated", "n_distinct_len"],
                                        ascending=[True, False])
    with pd.option_context("display.width", 200, "display.max_rows", 60):
        print()
        print("=== indel-length diversity at called sites ===")
        print(df.to_string(index=False))
    print()
    for g, sub in df.groupby("curated"):
        print(f"{g:9s} n={len(sub):2d}  median distinct lengths={sub['n_distinct_len'].median():.1f}"
              f"  max={sub['n_distinct_len'].max()}"
              f"  median modal_len_frac={sub['modal_len_frac'].median():.3f}")
    print()
    print("Position should be concordant; length should not. A site whose every read carries the")
    print("identical indel is more suspicious than one with varied lengths at the same cut.")

    if args.out:
        df.to_csv(args.out, index=False)
        print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
