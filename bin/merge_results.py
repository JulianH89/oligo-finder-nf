#!/usr/bin/env python
import sys
import argparse
import pandas as pd

def load_data(file_path):
    try:
        return pd.read_csv(file_path, sep="\t")
    except Exception as e:
        print(f"Error loading file {file_path}: {e}", file=sys.stderr)
        sys.exit(1)

def merge_results(filtered_metadata_path, crossreactivity_report_path, target_accessibility_path, output_path):
    # Load the filtered metadata and cross-reactivity report
    filtered_metadata = load_data(filtered_metadata_path)
    crossreactivity_report = load_data(crossreactivity_report_path)
    target_accessibility = load_data(target_accessibility_path)

    # Merge the two DataFrames on the 'ID' column
    merged = pd.merge(filtered_metadata, crossreactivity_report, on="#ID", how="right")
    
    # Merge with target accessibility data
    merged = pd.merge(merged, target_accessibility, on="#ID", how="left")

    # Save the merged DataFrame to the output file
    merged.to_csv(output_path, sep="\t", index=False)


def main():
    parser = argparse.ArgumentParser(description="Merge filtered metadata and cross-reactivity reports.")
    parser.add_argument("--filtered_metadata", required=True, help="Path to the filtered metadata TSV file.")
    parser.add_argument("--crossreactivity_report", required=True, help="Path to the cross-reactivity report TSV file.")
    parser.add_argument("--target_accessibility", help="Path to the target accessibility TSV file.")
    parser.add_argument("--output", required=True, help="Path to the output merged TSV file.")
    args = parser.parse_args()

    merge_results(args.filtered_metadata, args.crossreactivity_report, args.target_accessibility, args.output)
    

if __name__ == "__main__":
    main()