#!/usr/bin/env python3
"""Convert a combined off-target sites CSV into the per-base hotspot VCF the OFFTARGET
arm consumes as a ``target_file``.

Each predicted site (a ``Start`` on a ``Chromosome`` in the combined sites table) is
expanded to a ``±window`` bp interval; overlapping intervals are merged; and every base in
the merged intervals becomes a VCF record ``CHROM POS . <REF> N . PASS`` with the reference
base fetched from the FASTA. Records are written with a contig header and sorted to the
reference's contig order, so the result opens cleanly with ``pysam.VariantFile`` in
``bin/find_edited_reads.py``.

This is the productionized successor to ``workflow/make_hotspot_file/make_hotspot_file.py``
(hardcoded FASTA path replaced by ``--fasta``; robust pysam writer from
``bin/make_hotspot_vcf.py``).

    targets_csv_to_vcf.py --csv guide.targets.csv --fasta ref.fa --window 200 -o guide.targets.vcf
"""
import argparse
import os
import sys

import pandas as pd
import pysam

__version__ = "1.0.0"


def check_file(path):
    if not os.path.isfile(path):
        raise argparse.ArgumentTypeError(f"File not found: {path}")
    return path


def merged_intervals(df):
    """[(chrom, start, end)] merged per chrom; inputs are 0-based half-open [start, end)."""
    out = []
    for chrom, grp in df.sort_values(["Chromosome", "Start"]).groupby("Chromosome", sort=False):
        cur_s = cur_e = None
        for s, e in zip(grp["Start"], grp["End"]):
            if cur_s is None:
                cur_s, cur_e = s, e
            elif s <= cur_e:                       # overlap or touch -> extend
                cur_e = max(cur_e, e)
            else:
                out.append((chrom, cur_s, cur_e))
                cur_s, cur_e = s, e
        if cur_s is not None:
            out.append((chrom, cur_s, cur_e))
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--csv", type=check_file, help="Combined off-target sites CSV (needs Chromosome,Start).")
    ap.add_argument("--bed", type=check_file, help="Optional BED of extra hotspot regions.")
    ap.add_argument("--fasta", type=check_file, required=True, help="Indexed reference FASTA.")
    ap.add_argument("--window", type=int, default=200, help="bp up/downstream of each site.")
    ap.add_argument("-o", "--outfile", required=True, help="Output VCF.")
    ap.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    args = ap.parse_args(argv)

    if not (args.csv or args.bed):
        ap.error("provide --csv and/or --bed")

    regions = []
    if args.csv:
        csv = pd.read_csv(args.csv, usecols=["Chromosome", "Start"])
        csv = csv.dropna(subset=["Chromosome", "Start"])
        csv["Start"] = csv["Start"].astype(int)
        # window around the site; clamp lower bound at 0 (0-based half-open)
        r = pd.DataFrame({
            "Chromosome": csv["Chromosome"],
            "Start": (csv["Start"] - args.window).clip(lower=0),
            "End": csv["Start"] + args.window,
        })
        regions.append(r)
    if args.bed:
        bed = pd.read_csv(args.bed, sep="\t", usecols=[0, 1, 2],
                          names=["Chromosome", "Start", "End"])
        regions.append(bed)

    df = pd.concat(regions, ignore_index=True)
    if df.empty:
        print("WARNING: no target regions; writing header-only VCF", file=sys.stderr)

    fasta = pysam.FastaFile(args.fasta)
    valid = set(fasta.references)

    # Build the VCF header with contigs from the FASTA (needed for a sorted, indexable VCF).
    header = pysam.VariantHeader()
    for contig in fasta.references:
        header.contigs.add(contig, length=fasta.get_reference_length(contig))

    # Expand merged intervals to per-base records, collected then sorted to contig order.
    records = []  # (contig_rank, pos, chrom, ref_base)
    rank = {c: i for i, c in enumerate(fasta.references)}
    for chrom, start, end in merged_intervals(df):
        if chrom not in valid:
            print(f"WARNING: {chrom} not in reference; skipping", file=sys.stderr)
            continue
        seq = fasta.fetch(chrom, start, end)
        for i, base in enumerate(seq):
            records.append((rank[chrom], start + i, chrom, base))
    records.sort(key=lambda t: (t[0], t[1]))

    with pysam.VariantFile(args.outfile, "w", header=header) as vcf_out:
        for _, pos0, chrom, base in records:
            rec = vcf_out.new_record()
            rec.chrom = chrom
            rec.pos = pos0 + 1                      # VCF is 1-based; pos0 is 0-based
            rec.id = "."
            rec.ref = base if base and base.upper() in "ACGTN" else "N"
            rec.alts = ("N",)
            rec.filter.add("PASS")
            vcf_out.write(rec)
    return 0


if __name__ == "__main__":
    sys.exit(main())
