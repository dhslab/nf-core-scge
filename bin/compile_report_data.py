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
    parser.add_argument("--cna_plot", required=False, default=None, help="Path to CNA plot PNG.")
    parser.add_argument("--baf_plot", required=False, default=None, help="Path to BAF plot PNG.")
    parser.add_argument("--circos_plot", required=False, default=None, help="Path to Circos plot PNG.")
    parser.add_argument("--on_target_sv_transgene", required=True, help="Path to VEP-annotated on-target SV and transgene integration TSV.")
    parser.add_argument("--vcf_tsv", required=True, help="Path to VEP-annotated small variant TSV for targeted gene mutations.")
    parser.add_argument("--off_target_indels", required=True, help="Path to off-target indel analysis file.")
    parser.add_argument("--control_sample", required=False, default="N/A", help="Control/normal sample identifier.")
    parser.add_argument("--grnas", required=False, default="", help="Comma-separated list of gRNAs.")
    parser.add_argument("--coverage_metrics", action='append', default=None, help="Path(s) to coverage metrics files (can pass multiple).")
    parser.add_argument("-o", "--output", required=True, help="Output JSON file path.")
    
    args = parser.parse_args()

    # Parse the off-target indel file
    off_target_data = parse_offtarget_file(args.off_target_indels)

    # Parse the on-target SV and transgene data
    # Keep 'NA' as literal strings to avoid NaN in JSON, and read all columns as strings
    try:
        on_target_sv_transgene_df = pd.read_csv(
            args.on_target_sv_transgene,
            sep='	',
            comment='#',
            keep_default_na=False,
            na_filter=False,
            dtype=str,
        )
    except pd.errors.EmptyDataError:
        on_target_sv_transgene_df = pd.DataFrame(columns=[
            "Location", "Consequence", "SYMBOL", "BIOTYPE", "EXON",
            "INTRON", "STRAND", "Canonical", "Pick"
        ])
    on_target_sv_transgene_data = on_target_sv_transgene_df.to_dict(orient='records')

    # Parse targeted gene mutations (small variants) TSV if present
    targeted_gene_mutations = []
    try:
        sv_df = pd.read_csv(
            args.vcf_tsv,
            sep='	',
            comment='#',
            keep_default_na=False,
            na_filter=False,
            dtype=str,
        )
        # Summarize results for canonical list; fall back gracefully
        genes_of_interest = ["TP53", "DNMT3A", "RUNX1", "TET2"]
        for gene in genes_of_interest:
            has_variant = False
            # Try common columns
            gene_col = None
            for cand in ["SYMBOL", "Gene", "gene", "Symbol", "symbol"]:
                if cand in sv_df.columns:
                    gene_col = cand
                    break
            consequence_col = None
            for cand in ["Consequence", "CSQ", "consequence"]:
                if cand in sv_df.columns:
                    consequence_col = cand
                    break
            if gene_col is not None:
                rows = sv_df[sv_df[gene_col] == gene]
                if not rows.empty:
                    has_variant = True
                    if consequence_col is not None:
                        effects = sorted(set(
                            ",".join(rows[consequence_col].astype(str)).split(",")
                        ))
                        targeted_gene_mutations.append({"gene": gene, "result": "; ".join([e for e in effects if e])})
                    else:
                        targeted_gene_mutations.append({"gene": gene, "result": "variant detected"})
            if not has_variant:
                targeted_gene_mutations.append({"gene": gene, "result": "no mutations identified"})
    except Exception:
        # If TSV missing or unreadable, emit default rows
        for gene in ["TP53", "DNMT3A", "RUNX1", "TET2"]:
            targeted_gene_mutations.append({"gene": gene, "result": "no data"})

    # Create a dictionary to hold all the report data.
    report_data = {
        "sample_id": args.sample_id,
        "transgene_description": args.transgene,
        "plots": {
            "cna": args.cna_plot,
            "baf": args.baf_plot,
        },
        "tables": {
            "on_target_sv_transgene": on_target_sv_transgene_data,
            "off_target_indels": off_target_data,
            "targeted_gene_mutations": targeted_gene_mutations
        },
        "metadata": {
            "drug_product": args.sample_id,
            "control_sample": args.control_sample,
            "assay": "WGS",
            "grnas": [s.strip() for s in args.grnas.split(",") if s.strip()],
            "mean_coverage": { "tumor": None, "normal": None }
        }
    }

    if args.circos_plot:
        report_data["plots"]["circos"] = args.circos_plot

    # Parse coverage metrics if provided
    def try_parse_mean_coverage(file_path):
        try:
            # Heuristic: support simple TSV/CSV with columns including 'mean_coverage' or 'MEAN_COVERAGE'
            import os
            import csv
            with open(file_path, 'r') as fh:
                sample_content = fh.read(4096)
            delimiter = '\t' if '\t' in sample_content and ',' not in sample_content else ','
            rows = []
            with open(file_path, 'r') as fh:
                reader = csv.DictReader(fh, delimiter=delimiter)
                for row in reader:
                    rows.append({k.strip(): v for k, v in row.items()})
            if not rows:
                return None
            header = rows[0].keys()
            cov_key = None
            for cand in ['mean_coverage','MEAN_COVERAGE','Mean_Coverage','MEAN_COV','MEAN']:
                if cand in header:
                    cov_key = cand
                    break
            sample_key = None
            for cand in ['sample','SAMPLE','id','ID','name','NAME']:
                if cand in header:
                    sample_key = cand
                    break
            if cov_key is None:
                return None
            # If two rows exist, assume first is tumor (case id), second is control if present
            tumor_cov = None
            normal_cov = None
            if sample_key is not None:
                for r in rows:
                    sid = str(r[sample_key])
                    if sid == args.sample_id:
                        tumor_cov = float(r[cov_key])
                    elif sid == args.control_sample:
                        normal_cov = float(r[cov_key])
            if tumor_cov is None and rows:
                tumor_cov = float(rows[0][cov_key])
            if normal_cov is None and len(rows) > 1:
                normal_cov = float(rows[1][cov_key])
            return {"tumor": tumor_cov, "normal": normal_cov}
        except Exception:
            return None

    if args.coverage_metrics:
        for cov_path in args.coverage_metrics:
            cov = try_parse_mean_coverage(cov_path)
            if cov is not None:
                # Only set if values exist; prefer first successful parse
                if cov.get('tumor') is not None:
                    report_data['metadata']['mean_coverage']['tumor'] = cov['tumor']
                if cov.get('normal') is not None:
                    report_data['metadata']['mean_coverage']['normal'] = cov['normal']
                if report_data['metadata']['mean_coverage']['tumor'] is not None and report_data['metadata']['mean_coverage']['normal'] is not None:
                    break

    # Minimal schema validation
    def require(path, container):
        cur = container
        for key in path:
            if isinstance(cur, dict) and key in cur:
                cur = cur[key]
            else:
                raise ValueError(f"Missing required key in report_data: {'/'.join(path)}")
        return cur

    # Validate critical keys
    require(["sample_id"], report_data)
    if args.cna_plot is not None:
        require(["plots","cna"], report_data)
    if args.baf_plot is not None:
        require(["plots","baf"], report_data)
    require(["tables","on_target_sv_transgene"], report_data)
    require(["tables","off_target_indels"], report_data)
    require(["metadata","drug_product"], report_data)
    require(["metadata","control_sample"], report_data)

    # Write the compiled data to the output JSON file.
    with open(args.output, 'w') as f:
        json.dump(report_data, f, indent=4, allow_nan=False)

if __name__ == "__main__":
    main() 