#!/usr/bin/env python3
"""
bnd_from_indels_to_vcf.py -- turn the breakend calls in an *.offtarget_analysis.tsv
into a VCF, one record per junction, for the per-sample SCGE report.

This is the UNFILTERED breakend set: every junction the caller kept survives to the
VCF. The triaged view is a separate, cohort-level product (review_filter_bnd.py ->
bnd_review_queue.tsv); the two share a source column and nothing else.

The bnd_info payload is 11 pipe-separated fields per junction, joined by ';',
written by find_edited_reads.py:

    chrom|pos|chrom2|pos2|strands|ref|alt|distance|distance2|counts|control_alt_counts
      0    1     2     3     4     5   6      7         8        9          10

Fields 5 and 6 are already a valid VCF REF and BND ALT -- find_edited_reads.py's
format_bnd() builds them, brackets and orientation included (e.g. "CCAC]chr2:32916560]"
or "]chr2:32916405]GCC"). We use them verbatim rather than re-deriving the bracket
convention from `strands`, which is both simpler and keeps the VCF consistent with the
caller that produced it.

Mates are LINKED, not synthesized: bnd_info already carries a junction from both ends
as separate records, so the reciprocal partner is normally present in the same file and
carries its own correctly-oriented ALT. We pair them by coordinate and emit MATEID.
"""
import argparse
import csv
import sys

# A pysam-free reader: this runs in docker-baseimage, which is a thinner image than the
# review processes use.
csv.field_size_limit(1 << 30)   # bnd_info is one long field; the 128 KB default truncates it


def parse_junctions(path):
    """Yield (chrom, pos, chrom2, pos2, strands, ref, alt, counts, control) per junction."""
    seen = set()
    out = []
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            info = (row.get("bnd_info") or ".").strip()
            if not info or info == ".":
                continue
            for ev in info.split(";"):
                f = ev.split("|")
                if len(f) < 11:
                    continue
                try:
                    chrom, pos, chrom2, pos2 = f[0], int(f[1]), f[2], int(f[3])
                    counts, control = int(f[9]), int(f[10])
                except (ValueError, IndexError):
                    continue
                ref, alt = f[5] or "N", f[6]
                if not alt or alt == ".":
                    continue
                # The same junction is reported from each end, and adjacent target sites
                # can report it twice; collapse exact repeats so one junction is one record.
                key = (chrom, pos, chrom2, pos2, f[4], alt)
                if key in seen:
                    continue
                seen.add(key)
                out.append((chrom, pos, chrom2, pos2, f[4], ref, alt, counts, control))
    return out


def link_mates(records):
    """Map record index -> mate index, pairing a junction with its reciprocal record."""
    by_end = {}
    for i, r in enumerate(records):
        by_end.setdefault((r[0], r[1], r[2], r[3]), i)
    mate = {}
    for i, r in enumerate(records):
        j = by_end.get((r[2], r[3], r[0], r[1]))
        if j is not None and j != i:
            mate[i] = j
    return mate


HEADERS = [
    "##fileformat=VCFv4.2",
    '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
    '##INFO=<ID=MATEID,Number=1,Type=String,Description="ID of mate breakend">',
    '##INFO=<ID=EVENT,Number=1,Type=String,Description="ID of associated event">',
    '##INFO=<ID=SR,Number=1,Type=Integer,Description="Reads supporting this junction in the edited sample">',
    '##INFO=<ID=CTRL,Number=1,Type=Integer,Description="Reads supporting this junction in the matched control">',
]


def main(args):
    records = parse_junctions(args.indels_path)
    mate = link_mates(records)

    # Coordinate-sorted output, but IDs are assigned before sorting so a mate's ID does
    # not depend on sort order.
    ids = {i: f"BND{i + 1}" for i in range(len(records))}
    order = sorted(range(len(records)), key=lambda i: (records[i][0], records[i][1]))

    with open(args.outfile, "w") as out:
        for h in HEADERS:
            out.write(h + "\n")
        out.write(f"##sampleId={args.meta_id}\n")
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for i in order:
            chrom, pos, chrom2, pos2, _strands, ref, alt, counts, control = records[i]
            info = [f"SVTYPE=BND", f"SR={counts}", f"CTRL={control}"]
            if i in mate:
                info.append(f"MATEID={ids[mate[i]]}")
                info.append(f"EVENT=EVT{min(i, mate[i]) + 1}")
            out.write(f"{chrom}\t{pos}\t{ids[i]}\t{ref}\t{alt}\t.\tPASS\t{';'.join(info)}\n")

    print(f"{args.meta_id}: {len(records)} breakend records -> {args.outfile}", file=sys.stderr)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert an offtarget analysis table to a BND VCF.")
    parser.add_argument("--meta_id", required=True, help="Sample ID")
    parser.add_argument("--indels_path", required=True, help="Path to *.offtarget_analysis.tsv")
    parser.add_argument("--outfile", required=True, help="Path to output VCF file")
    main(parser.parse_args())
