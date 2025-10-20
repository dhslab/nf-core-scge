#!/usr/bin/env python3

import argparse
import sys
import pandas as pd
import pyranges as pr
import pysam
from pathlib import Path
from typing import List

# Define the script version
__version__ = "1.0.0"

# --- Helper Functions ---

def check_file(path: str) -> Path:
    """Checks if a file path exists and is a file."""
    p = Path(path)
    if not p.is_file():
        raise argparse.ArgumentTypeError(f"File not found: {path}")
    return p

def add_sequence_column(row: pd.Series, fasta_handle: pysam.FastaFile) -> str:
    """
    Fetches a genomic sequence for a given row using an open pysam.FastaFile handle.
    """
    try:
        # pysam.fetch is 0-based, and BED Start is 0-based.
        # It fetches the interval [start, end), which matches the BED standard.
        sequence = fasta_handle.fetch(row["Chromosome"], row["Start"], row["End"])
        return sequence
    except (ValueError, KeyError) as e:
        print(f"Warning: Could not fetch sequence for {row['Chromosome']}:{row['Start']}-{row['End']}. Reason: {e}", file=sys.stderr)
        return "N" * (row['End'] - row['Start'])


def dataframe_to_vcf(df: pd.DataFrame, fasta_handle: pysam.FastaFile) -> str:
    """
    Converts a DataFrame with genomic positions into a VCF formatted string,
    including a proper header with contigs from the reference FASTA.
    """
    vcf_lines = ["##fileformat=VCFv4.2"]

    # Add contig lines to the header from the FASTA index
    for contig in fasta_handle.references:
        length = fasta_handle.get_reference_length(contig)
        vcf_lines.append(f"##contig=<ID={contig},length={length}>")
        
    vcf_lines.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
    
    # Prepare DataFrame for VCF format
    df_vcf = df.copy()
    df_vcf['#CHROM'] = df_vcf['Chromosome']
    df_vcf['POS'] = df_vcf['Position'] + 1 # VCF is 1-based
    df_vcf['ID'] = '.'
    df_vcf['REF'] = df_vcf['Sequence']
    df_vcf['ALT'] = '.'  # No alternative allele for a hotspot VCF
    df_vcf['QUAL'] = '.'
    df_vcf['FILTER'] = 'PASS'
    df_vcf['INFO'] = '.'
    
    vcf_body = df_vcf[['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']].to_csv(sep='\t', index=False, header=False)
    
    return "\n".join(vcf_lines) + "\n" + vcf_body

# --- Main Application Logic ---

def main():
    """
    Main function to generate a hotspot VCF from a BED file of genomic regions.
    """
    parser = argparse.ArgumentParser(
        description="Prepare a hotspot VCF from a BED file of genomic regions.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--bed", type=check_file, required=True, help="BED file containing hotspot regions.")
    parser.add_argument("--fasta", type=check_file, required=True, help="Path to the indexed reference genome FASTA file.")
    parser.add_argument("--outfile", required=True, help="Output VCF file.")
    parser.add_argument('--version', action='version', version=f'%(prog)s {__version__}')
    args = parser.parse_args()

    try:
        # --- 1. Open FASTA file once ---
        print(f"Opening reference FASTA: {args.fasta}", file=sys.stderr)
        fasta_handle = pysam.FastaFile(str(args.fasta))

        # --- 2. Read and Merge Genomic Regions ---
        print(f"Reading and processing BED file: {args.bed}", file=sys.stderr)
        bed_df = pd.read_csv(args.bed, sep='\t', usecols=[0, 1, 2], names=['Chromosome', 'Start', 'End'])
        merged_df = pr.PyRanges(bed_df).merge().sort().df

        # --- 3. Fetch Sequences for Merged Regions ---
        print("Fetching sequences for merged regions...", file=sys.stderr)
        merged_df['sequences'] = merged_df.apply(
            lambda row: add_sequence_column(row, fasta_handle), 
            axis=1
        )

        # --- 4. Expand Regions into a VCF-like DataFrame (Optimized) ---
        print("Expanding genomic regions into VCF format...", file=sys.stderr)
        all_chroms, all_positions, all_sequences = [], [], []

        for _, row in merged_df.iterrows():
            num_bases = len(row['sequences'])
            all_chroms.extend([row['Chromosome']] * num_bases)
            all_positions.extend(range(row['Start'], row['End']))
            all_sequences.extend(list(row['sequences']))

        vcf_df = pd.DataFrame({
            'Chromosome': all_chroms,
            'Position': all_positions,
            'Sequence': all_sequences
        })

        # --- 5. Convert to VCF and Write to File ---
        result_vcf_str = dataframe_to_vcf(vcf_df, fasta_handle)
        output_filename = f"{args.outfile}"
    
        with open(output_filename, "w") as f:
            f.write(result_vcf_str)
        print(f"Successfully wrote hotspot VCF to '{output_filename}'", file=sys.stderr)

    except Exception as e:
        sys.exit(f"An error occurred: {e}")
    finally:
        # Ensure the FASTA file handle is closed
        if 'fasta_handle' in locals() and fasta_handle:
            fasta_handle.close()

if __name__ == "__main__":
    main()

