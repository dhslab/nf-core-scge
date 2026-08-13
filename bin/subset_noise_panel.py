#!/usr/bin/env python3
"""
subset_noise_panel.py — cut a DRAGEN systematic-noise panel down to the sites a guide can query.

Why this exists
---------------
The IDPF SNV panel is 1 GB compressed / 97,851,753 records, and `load_snv_noise` streams it
whole because it ships without a tabix index. Paid once per cohort that is ~3 minutes and nobody
notices. Paid once per SAMPLE -- which is what happens if the panel lookup moves into the caller --
it becomes ~3 minutes x N samples, and on the 32-sample CAR-T cohort that is over an hour of pure
re-reading of the same file.

The fix is that the panel only ever gets asked about positions in the guide's target file, and a
target file is per-GUIDE, not per-sample. So subset once per guide panel and every sample sharing
that guide reads a file that fits in memory.

Measured on the 32-sample CAR-T cohort: 97,851,753 records -> 4,891 (60 KB), a 20,000x reduction,
and `review_filter.py` produces byte-identical output from the subset. The reduction is lossless by
construction: `load_snv_noise` keeps a record only when its END column lands within +/-slop of a
queried position, so anything dropped here could never have been returned.

This also retires the tabix question -- an index solves a problem the subset removes outright.

usage:
  # from a caller output table (uses its chrom/end columns)
  subset_noise_panel.py --panel IDPF_...snv.bed.gz --sites SAMPLE.offtarget_analysis.tsv \\
                        -o panel.subset.bed.gz

  # from the guide's target file, before any sample has been called
  subset_noise_panel.py --panel IDPF_...snv.bed.gz --targets GUIDE.targets.csv \\
                        -o panel.subset.bed.gz
"""
import argparse
import glob
import gzip
import sys

import pandas as pd

DEFAULT_SLOP = 2


def positions_from_tables(paths, slop):
    """Union of (chrom, end +/- slop) over caller output tables."""
    want = {}
    n = 0
    for p in paths:
        df = pd.read_csv(p, sep="\t", usecols=["chrom", "end"])
        n += len(df)
        for c, e in zip(df.chrom, df.end):
            s = want.setdefault(str(c), set())
            for d in range(-slop, slop + 1):
                s.add(int(e) + d)
    return want, n


def positions_from_targets(paths, slop, window):
    """Union of target windows from a guide target file (CSV/BED/VCF-ish).

    A target file gives the site, not the exact indel position, so each target contributes a
    window rather than a point -- `window` should be at least the caller's --target-window.
    """
    want = {}
    n = 0
    for p in paths:
        sep = "\t" if p.endswith((".bed", ".tsv", ".txt")) else ","
        df = pd.read_csv(p, sep=sep, comment="#")
        cols = {c.lower(): c for c in df.columns}
        cc = cols.get("chrom") or cols.get("chromosome") or cols.get("#chrom") or df.columns[0]
        pc = cols.get("start") or cols.get("pos") or cols.get("position") or df.columns[1]
        n += len(df)
        for c, p0 in zip(df[cc], df[pc]):
            try:
                p0 = int(p0)
            except (TypeError, ValueError):
                continue
            s = want.setdefault(str(c), set())
            for d in range(-window - slop, window + slop + 1):
                s.add(p0 + d)
    return want, n


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--panel", required=True, help="systematic_noise.snv.bed.gz")
    ap.add_argument("--sites", nargs="*", default=[], help="*.offtarget_analysis.tsv")
    ap.add_argument("--targets", nargs="*", default=[], help="guide target file(s)")
    ap.add_argument("-o", "--out", required=True, help="output .bed.gz")
    ap.add_argument("--slop", type=int, default=DEFAULT_SLOP,
                    help=f"must be >= the slop the consumer uses (default {DEFAULT_SLOP}); too "
                         f"small silently drops records the filter would have matched")
    ap.add_argument("--window", type=int, default=150,
                    help="half-width around each target site for --targets (default 150, the "
                         "caller's --target-window)")
    a = ap.parse_args()

    paths = lambda pats: [f for p in pats
                          for f in (sorted(glob.glob(p)) if any(ch in p for ch in "*?[") else [p])]
    if a.sites:
        want, n_in = positions_from_tables(paths(a.sites), a.slop)
        src = f"{n_in} table rows"
    elif a.targets:
        want, n_in = positions_from_targets(paths(a.targets), a.slop, a.window)
        src = f"{n_in} target sites (+/-{a.window} bp)"
    else:
        sys.exit("ERROR: give --sites or --targets")

    n_pos = sum(len(v) for v in want.values())
    print(f"queryable positions: {n_pos} from {src}")

    kept = seen = 0
    op = gzip.open if a.panel.endswith(".gz") else open
    with op(a.panel, "rt") as fh, gzip.open(a.out, "wt") as og:
        for line in fh:
            if line[0] == "#":
                og.write(line)          # keep the header: ##PON SAMPLES gives N for --panel-p
                continue
            seen += 1
            f = line.split("\t", 3)
            s = want.get(f[0])
            if s is None:
                continue
            try:
                if int(f[2]) in s:
                    og.write(line)
                    kept += 1
            except (ValueError, IndexError):
                continue
    pct = kept / seen * 100 if seen else 0
    print(f"panel records {seen} -> {kept} ({pct:.4f}%)  -> {a.out}")


if __name__ == "__main__":
    main()
