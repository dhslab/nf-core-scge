#!/usr/bin/env python

# compile_report_data.py
# This script compiles data from various sources into a single JSON file for the Quarto report.

import argparse
import json
import pandas as pd

def parse_offtarget_file(file_path):
    """Parse the off-target file and return a list of dictionaries."""
    records = []
    with open(file_path) as data:
        header = data.readline().strip().split('\t')
        for line in data:
            stripped_line = line.strip().split("\t")
            if len(stripped_line) < len(header):
                continue
            
            row_dict = dict(zip(header, stripped_line))
            records.append(row_dict)
    return records

def main():
    """
    Main function to parse arguments and compile data.
    """
    parser = argparse.ArgumentParser(description="Compile data for Quarto report.")
    parser.add_argument("--sample_id", required=True, help="Sample ID.")
    parser.add_argument("--transgene", required=True, help="Transgene description string.")
    parser.add_argument("--cna_plot", required=True, help="Path to CNA plot PNG.")
    parser.add_argument("--baf_plot", required=True, help="Path to BAF plot PNG.")
    parser.add_argument("--circos_plot", required=True, help="Path to Circos plot PNG.")
    parser.add_argument("--on_target_sv_transgene", required=True, help="Path to VEP-annotated on-target SV and transgene integration TSV.")
    parser.add_argument("--off_target_indels", required=True, help="Path to off-target indel analysis file.")
    parser.add_argument("-o", "--output", required=True, help="Output JSON file path.")
    
    args = parser.parse_args()

    # Parse the off-target indel file
    off_target_data = parse_offtarget_file(args.off_target_indels)

    # Parse the on-target SV and transgene data
    on_target_sv_transgene_df = pd.read_csv(args.on_target_sv_transgene, sep='\t')
    on_target_sv_transgene_data = on_target_sv_transgene_df.to_dict(orient='records')


    # Create a dictionary to hold all the report data.
    report_data = {
        "sample_id": args.sample_id,
        "transgene_description": args.transgene,
        "plots": {
            "cna": args.cna_plot,
            "baf": args.baf_plot,
            "circos": args.circos_plot
        },
        "tables": {
            "on_target_sv_transgene": on_target_sv_transgene_data,
            "off_target_indels": off_target_data
        },
        "metadata": {
            # We can add other metadata here if needed.
        }
    }

    # Write the compiled data to the output JSON file.
    with open(args.output, 'w') as f:
        json.dump(report_data, f, indent=4)

if __name__ == "__main__":
    main() 