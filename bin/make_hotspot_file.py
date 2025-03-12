#!/usr/bin/env python3

import sys, pysam, argparse, os
import pandas as pd, pyranges as pr

def checkfile(file_path):
    """Check if a file exists at the given path."""
    if not os.path.exists(file_path):
        raise argparse.ArgumentTypeError(f"The file {file_path} does not exist.")
    return file_path

def get_position(fasta_file, chromosome, start, end):
    fasta = pysam.FastaFile(fasta_file)
    sequences = fasta.fetch(chromosome, start, end)
    return sequences

def dataframe_to_vcf(df):
    vcf_lines = []
    vcf_lines.append("##fileformat=VCFv4.1")
    vcf_lines.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
    
    for _, row in df.iterrows():
        chromosome = row["Chromosome"]
        position = row["Position"]
        sequence = row["Sequence"]
        
        vcf_line = f"{chromosome}\t{position}\t.\t{sequence}\tN\t.\t.\t."
        vcf_lines.append(vcf_line)
    
    return "\n".join(vcf_lines)

def add_sequence_column(row):
    chromosome = row["Chromosome"]
    start = row["Start"]
    end = row["End"]
    fasta_file = "/storage1/fs1/dspencer/Active/spencerlab/refdata/hg38/all_sequences.fa"
    sequences = get_position(fasta_file, chromosome, start, end+1)
    return sequences

def main():
    parser = argparse.ArgumentParser(description="Prepare hotspot vcf from sites file and/or bed file")
    parser.add_argument('--bed', type=checkfile, help="bed hotspot file")
    parser.add_argument('--csv', type=checkfile, help="csv hotspot file")
    parser.add_argument('--window', type=int, default=200, help="window")
    parser.add_argument('--id', help="meta id")
    args = parser.parse_args()

    window = args.window 
    data_df = pd.DataFrame()
    if (args.bed):
        bed_file = os.path.realpath(args.bed)
        bed_df = pd.read_csv(bed_file, sep='\t', usecols=[0, 1, 2], names=['Chromosome', 'Start', 'End'])
        data_df = pd.concat([data_df, bed_df], ignore_index=True)
    if (args.csv):
        csv_file = os.path.realpath(args.csv)
        csv_df = pd.read_csv(csv_file, usecols=["Chromosome", "Start"])
        csv_df['End'] = csv_df['Start'] + window
        csv_df['Start'] = csv_df['Start'].apply(lambda x: x - window)
        data_df = pd.concat([data_df, csv_df], ignore_index=True)

    data_p = pr.PyRanges(data_df).merge().sort()
    df = data_p.df
    df['sequences'] = df.apply(add_sequence_column, axis=1)
    
    vcf_df = pd.DataFrame()
    for _,rows in df.iterrows():
        data = {'Position': list(range(rows['Start'],rows['End']+1)), 'Sequence': list(rows['sequences']) }
        row_df = pd.DataFrame(data)
        row_df['Chromosome'] = rows['Chromosome']
        vcf_df = pd.concat([vcf_df, row_df], ignore_index=True)
    result_vcf = dataframe_to_vcf(vcf_df)
    id = args.id
    with open(f"{id}.hotspot.vcf", "w") as f:
        f.write(result_vcf)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python script.py --bed <bed_file> --csv <csv_file> --window [window] --id <id>")
        sys.exit(1)
    else:
        main()
