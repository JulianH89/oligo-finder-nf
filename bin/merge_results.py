#!/usr/bin/env python
import sys
import argparse
import pandas as pd

def load_data(file_path):
    try:
        data = pd.read_csv(file_path, sep="\t")
        return data
    except Exception as e:
        print(f"Error loading file {file_path}: {e}", file=sys.stderr)
        sys.exit(1)

def merge_results(filtered_metadata_path, crossreactivity_report_path, output_path):
    # Load the filtered metadata and cross-reactivity report
    filtered_metadata = load_data(filtered_metadata_path)
    crossreactivity_report = load_data(crossreactivity_report_path)

    # Merge the two DataFrames on the 'ID' column
    merged = pd.merge(filtered_metadata, crossreactivity_report, on="#ID", how="right")

    # Save the merged DataFrame to the output file
    merged.to_csv(output_path, sep="\t", index=False)


def main():
    parser = argparse.ArgumentParser(description="Merge filtered metadata and cross-reactivity reports.")
    parser.add_argument("--filtered_metadata", required=True, help="Path to the filtered metadata TSV file.")
    parser.add_argument("--crossreactivity_report", required=True, help="Path to the cross-reactivity report TSV file.")
    parser.add_argument("--output", required=True, help="Path to the output merged TSV file.")
    args = parser.parse_args()

    merge_results(args.filtered_metadata, args.crossreactivity_report, args.output)
    

if __name__ == "__main__":
    main()