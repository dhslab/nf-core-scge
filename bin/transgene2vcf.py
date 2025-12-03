#!/usr/bin/env python3

import sys
import csv
import argparse
import os

__version__ = "1.0.0"

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Convert a transgene TSV file into a VCF file."
    )
    # Positional arguments do not use -- dashes
    parser.add_argument(
        "transgene_path", 
        help="Path to the input transgene TSV file."
    )
    parser.add_argument(
        "output_file", 
        help="Path for the output VCF file."
    )
    parser.add_argument(
        "-v", "--version",
        action="version",
        version=f"%(prog)s {__version__}",
        help="Show the version number and exit."
    )
    return parser.parse_args()

def main():
    args = parse_args()

    transgene_path = args.transgene_path
    vcf_path = args.output_file

    # Check if input file exists
    if not os.path.exists(transgene_path):
        sys.exit(f"Error: Input file '{transgene_path}' not found.")

    records = []

    try:
        with open(transgene_path, 'r') as fh:
            reader = csv.reader(fh, delimiter='\t')
            
            # Read header safely
            try:
                header = next(reader)
            except StopIteration:
                sys.exit(f"Error: Input file '{transgene_path}' is empty.")

            # Create a map of column name to index
            hdr = {k: i for i, k in enumerate(header)}

            # Validate required columns exist
            required_cols = ['Chromosome', 'Start']
            for col in required_cols:
                if col not in hdr:
                    sys.exit(f"Error: Required column '{col}' missing from input file header.")

            for line_num, line in enumerate(reader, start=2):
                if not line: continue # Skip empty lines
                
                try:
                    chrom = line[hdr['Chromosome']]
                    pos = line[hdr['Start']]
                    ref = 'N'
                    alt = '<INS>'
                    records.append((chrom, pos, ref, alt))
                except IndexError:
                    print(f"Warning: Skipping malformed line {line_num} in {transgene_path}", file=sys.stderr)

    except Exception as e:
        sys.exit(f"Error reading input file: {e}")

    # Write VCF
    try:
        with open(vcf_path, 'w') as out:
            out.write("##fileformat=VCFv4.2\n")
            out.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n')
            out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

            for idx, (chrom, pos, ref, alt) in enumerate(records, start=1):
                variant_id = f"TRANSGENE{idx}"
                qual = "."
                filter_val = "PASS"
                info = "SVTYPE=TRANSGENE"
                out.write(f"{chrom}\t{pos}\t{variant_id}\t{ref}\t{alt}\t{qual}\t{filter_val}\t{info}\n")
        
    except IOError as e:
        sys.exit(f"Error writing to output file '{vcf_path}': {e}")

if __name__ == '__main__':
    main()