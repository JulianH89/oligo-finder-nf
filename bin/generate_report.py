#!/usr/bin/env python

import argparse
import json
import csv
import sys

def generate_report(json_file_path, output_tsv_path):
    """
    Reads a JSON file from the parse_sam step and generates a final
    tab-separated report with the specified format.
    """
    try:
        with open(json_file_path, 'r') as f:
            parsed_data = json.load(f)
    except (FileNotFoundError, json.JSONDecodeError) as e:
        print(f"Error reading or parsing JSON file {json_file_path}: {e}", file=sys.stderr)
        sys.exit(1)

    with open(output_tsv_path, 'w', newline='') as f_out:
        writer = csv.writer(f_out, delimiter='\t')

        # Write the new, corrected header
        writer.writerow(['Oligo_id', 'sequence', 'gc_content', 'mismatches', 'matched_accession', 'num_of_matched'])

        # Iterate through each oligo in the JSON data
        for oligo_id, data in parsed_data.items():
            sequence = data.get('sequence', 'N/A')
            gc_content = data.get('gc_content', 0.0)

            # Iterate through each mismatch level for the current oligo
            # The keys are strings in JSON, so we sort them numerically
            sorted_mismatch_levels = sorted(data.get('mismatches', {}).keys(), key=int)

            for mismatch_level in sorted_mismatch_levels:
                mismatch_info = data['mismatches'][mismatch_level]
                accessions = mismatch_info.get('accessions', [])
                num_of_matched = len(accessions)

                # Apply the rule for formatting the accession list
                if num_of_matched > 10:
                    matched_accession = 'too many to record'
                else:
                    matched_accession = ','.join(accessions)

                # Write the final row with the new column order
                writer.writerow([
                    oligo_id,
                    sequence,
                    gc_content,
                    mismatch_level, # This is the new 'mismatches' column
                    matched_accession,
                    num_of_matched
                ])

def main():
    parser = argparse.ArgumentParser(description="Generate a TSV report from a parsed SAM JSON file.")
    parser.add_argument("--json", required=True, help="Input JSON file path.")
    parser.add_argument("--output", required=True, help="Output TSV file path.")
    args = parser.parse_args()

    generate_report(args.json, args.output)

if __name__ == "__main__":
    main()