#!/usr/bin/env python3
"""Combine predicted off-target sites for one gRNA across tools into a single table.

Merges the native outputs of Cas-OFFinder (bulge mode), CRISPRme, and (optionally) IDT
into the canonical off-target sites schema the OFFTARGET arm consumes as a per-guide
``target_file``:

    Source,DNA_Sequence,PAM,Chromosome,Strand,Start,Bulge_Type,Mismatch,Bulge_Size,On_target

Sites predicted by more than one tool (same Chromosome+Start+PAM) collapse to one row
whose ``Source`` is the ``|``-joined list of contributing tools. This is the hardened,
argparse-driven successor to the ``Combine_offtarget_results.py`` prototype; at least one
source is required and any source may be omitted (IDT is off by default).

The per-source coordinate adjustments (strand-dependent, and different for each tool) are
load-bearing and preserved verbatim from the validated prototype.

    combine_offtarget_results.py \\
        --casoffinder guide.casoffinder.txt \\
        --crisprme    guide.crisprme.tsv \\
        [--idt        guide_IDT_off-target.xlsx] \\
        -o guide.targets.csv
"""
import argparse
import sys

import pandas as pd

__version__ = "1.0.0"

# The canonical schema (underscore headers, `Strand`) — matches assets/stub/*.targets.csv
# and what bin/make_hotspot_vcf.py / bin/find_edited_reads.py expect downstream.
COLUMNS = ["Source", "DNA_Sequence", "PAM", "Chromosome", "Strand", "Start",
           "Bulge_Type", "Mismatch", "Bulge_Size", "On_target"]


def parse_idt(path):
    """IDT off-target Excel export -> rows in the canonical schema.

    Expected columns: 'Sequence', 'PAM', '#MM', 'Locus' (e.g. 'chr19:+7900116').
    The on-target row has a blank/NaN '#MM'.
    """
    rows = []
    idt = pd.read_excel(path)
    for _, row in idt.iterrows():
        seq = str(row["Sequence"])
        pam = row["PAM"]
        mismatch = str(row["#MM"])              # 'nan' for the on-target row
        locus = str(row["Locus"]).split(":")
        chrom = locus[0]
        strand = locus[1][0]
        loci = int(locus[1][1:])
        # IDT reports the PAM-distal end on '+' and needs +1 on '-' to reach the site start
        adjusted = str(loci + len(seq)) if strand == "+" else str(loci + 1)
        on_target = 1 if mismatch == "nan" else 0
        rows.append(["IDT", seq, pam, chrom, strand, adjusted,
                     "mismatch", mismatch, "NA", on_target])
    return rows


def parse_casoffinder(path):
    """Cas-OFFinder v3 native bulge output (tab-sep, '#' header lines) -> canonical rows.

    Columns (v3): Id, Bulge Type, crRNA, DNA, Chromosome, Location, Direction, Mismatches, Bulge Size.
    Cas-OFFinder >=3.0 does DNA/RNA bulges natively (no separate wrapper) and prepends an `Id`
    column vs the old 2.4 / cas-offinder-bulge layout; its header/comment lines start with '#'.
    """
    rows = []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            bulge_type = c[1]
            dna = c[3]
            dna_only = dna[0:len(dna) - 4]
            pam = dna[len(dna) - 3:len(dna)]
            chrom = c[4]
            position = int(c[5])
            direction = c[6]
            adjusted = str(position + (len(dna) - 3)) if direction == "+" else str(position + 4)
            mismatches = c[7]
            bulge_size = c[8]
            on_target = 1 if (bulge_type == "X" and mismatches == "0" and bulge_size == "0") else 0
            rows.append(["CasOffFinder", dna_only, pam, chrom, direction, adjusted,
                         bulge_type, mismatches, bulge_size, on_target])
    return rows


def parse_crisprme(path):
    """CRISPRme (>=2.1) *_integrated_results.tsv (one header line) -> canonical rows.

    complete-search writes this loose in Results/<output>/. Columns used: Spacer+PAM[0],
    Chromosome[1], Start_coordinate[2], Strand[3], Aligned_protospacer+PAM_REF[5], PAM[7],
    Mismatches[8], Bulges[9], Bulge_type[15]. Cas-OFFinder and CRISPRme report a site's
    coordinate up to ~2 bp apart (bulge handling differs), so the same physical site may not
    merge to one row across tools — harmless, since the downstream targets VCF windows each
    site by +/- hotspot_window_size and merges overlapping intervals.
    """
    rows = []
    with open(path) as fh:
        fh.readline()                           # header (single line, no leading '#')
        for line in fh:
            if not line.strip():
                continue
            c = line.rstrip("\n").split("\t")
            bulge_type = c[15]
            dna = c[5]
            dna_only = dna[0:len(dna) - 3]
            pam = c[7]
            chrom = c[1]
            start = int(c[2])
            direction = c[3]
            adjusted = str(start + (len(dna) - 4)) if direction == "+" else str(start + 4)
            mismatches = c[8]
            bulge_size = c[9]
            on_target = 1 if (bulge_type == "X" and mismatches == "0" and bulge_size == "0") else 0
            rows.append(["CrisprME", dna_only, pam, chrom, direction, adjusted,
                         bulge_type, mismatches, bulge_size, on_target])
    return rows


def combine(rows):
    """Rows -> deduplicated DataFrame in canonical schema/column order."""
    df = pd.DataFrame(rows, columns=COLUMNS)

    # Numeric coercion so `min` is well-defined across mixed string inputs (IDT emits
    # 'nan'/'NA'); NaNs are ignored by min, and an all-NaN group stays NaN.
    df["Mismatch"] = pd.to_numeric(df["Mismatch"], errors="coerce")
    df["Bulge_Size"] = pd.to_numeric(df["Bulge_Size"], errors="coerce")

    # Collapse the same physical site predicted by multiple tools into one row.
    grouped = df.groupby(["Chromosome", "Start", "PAM"], as_index=False).agg({
        "Source": lambda x: "|".join(x),
        "DNA_Sequence": "first",
        "Strand": "first",
        "Bulge_Type": "first",
        "Mismatch": "min",
        "Bulge_Size": "min",
        "On_target": "max",
    })
    # Nullable integer display: whole numbers (0, 3, …) not floats (0.0), and an empty
    # cell for a genuinely-absent value (e.g. a pure-IDT on-target with no mismatch count).
    for col in ("Mismatch", "Bulge_Size"):
        grouped[col] = grouped[col].astype("Int64")
    return grouped[COLUMNS]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--casoffinder", help="Cas-OFFinder bulge output (.txt)")
    ap.add_argument("--crisprme", help="CRISPRme targets/best-hits (.tsv)")
    ap.add_argument("--idt", help="IDT off-target export (.xlsx); optional")
    ap.add_argument("--pam", help="PAM used for the search (recorded for provenance only)")
    ap.add_argument("-o", "--output", required=True, help="Output combined sites CSV")
    ap.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    args = ap.parse_args(argv)

    if not (args.casoffinder or args.crisprme or args.idt):
        ap.error("provide at least one of --casoffinder / --crisprme / --idt")

    rows = []
    if args.idt:
        rows += parse_idt(args.idt)
    if args.casoffinder:
        rows += parse_casoffinder(args.casoffinder)
    if args.crisprme:
        rows += parse_crisprme(args.crisprme)

    if not rows:
        print("WARNING: no off-target sites parsed from the provided source(s); "
              "writing header-only output", file=sys.stderr)
        pd.DataFrame(columns=COLUMNS).to_csv(args.output, index=False)
        return 0

    combine(rows).to_csv(args.output, index=False)
    return 0


if __name__ == "__main__":
    sys.exit(main())
