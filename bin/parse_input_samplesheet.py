#!/usr/bin/env python3

import argparse
import os

import pandas as pd

__version__ = "1.0.0"


def parse_arguments() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Process spreadsheet file for samples that require demultiplexing, alignment, or summary analysis."
    )
    parser.add_argument(
        "--input_file", "-i", required=True, help="Path to the input CSV file"
    )
    parser.add_argument(
        "--output_dir", "-o", help="Directory to save the output CSV files"
    )
    return parser.parse_args()


def clean_sample_name(value: str) -> str:
    """Remove whitespace and replace spaces with underscores."""
    if isinstance(value, str):
        return value.strip().replace(" ", "_")
    return value


def save_df(dataframe: pd.DataFrame, output_dir: str, filename: str) -> None:
    """Save dataframe to a CSV file, applying a cleaning function if provided."""
    if not dataframe.empty:
        dataframe.loc[:, "id"] = dataframe["id"].apply(
            clean_sample_name
        )

        dataframe = dataframe.dropna(axis=1, how="all")
        dataframe = dataframe.loc[:, dataframe.ne("").any()]

        dataframe.to_csv(f"{output_dir}/{filename}", index=False)


def process_input(input_file: str, output_dir: str) -> None:
    """Process input file, merge with sample info, and create output files."""
    ext = input_file.split(".")[-1]
    df = pd.read_csv(input_file, sep="\t" if ext == "tsv" else ",",engine="python",on_bad_lines="error")

    id_col = "id"
    if id_col not in df.columns:
        raise ValueError("No 'id' or 'Content_Desc' column found in input file!")

    # remove columns if not valid header
    valid_headers = [
        "id",
        "edited_id",
        "control_id",
        "edited_cram",
        "control_cram",
        "edited_bam",
        "control_bam",
        "edited_read1",
        "control_read1",
        "edited_read2",
        "control_read2",
        "fastq_list",
        "mgi_samplemap",
        "dragen_path"
    ]
    df = df.loc[:, [col for col in df.columns if col in valid_headers]]
    
    # check for duplicate id columns in df
    if df[id_col].duplicated().any():
        raise ValueError("Duplicate 'id' column found in input file!")
    
    # Summary analysis samples. This supercedes all other inputs--samples with dragen_path will only be analyzed.
    if "dragen_path" in df.columns:
        dragen_df = df.dropna(subset=["dragen_path"]).dropna(axis=1, how="all")
        save_df(dragen_df, output_dir, "analysis_samples.csv")
        df = df[~df.index.isin(dragen_df.index)] if not dragen_df.empty else df

    # Alignment samples. Note this includes samples to demux as well as already demuxed fastqs and cram/bam for realignment.
    alignment_cols = ["edited_cram",
        "control_cram",
        "edited_bam",
        "control_bam",
        "edited_read1",
        "control_read1",
        "edited_read2",
        "control_read2",
        "fastq_list",
        "mgi_samplemap"]
    
    alignment_df = df.dropna(
        subset=[col for col in alignment_cols if col in df], thresh=1
    )
    save_df(alignment_df, output_dir, "alignment_samples.csv")


def main() -> None:
    args = parse_arguments()

    output_dir = args.output_dir
    if not output_dir or not os.path.exists(output_dir):
        output_dir = os.getcwd()

    process_input(args.input_file, output_dir)


if __name__ == "__main__":
    main()
