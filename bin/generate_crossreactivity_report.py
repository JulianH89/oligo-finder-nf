#!/usr/bin/env python

import argparse
import json
import csv
import sys

def load_geneid_accession_map(geneid_accession_path):
    # sourcery skip: avoid-builtin-shadow
    """
    Load the GeneID to accession mapping from a given file.
    The file is expected to have two columns: GeneID and Accession, separated by tabs.
    """
    mapping = {}
    try:
        with open(geneid_accession_path, 'r') as f:
            next(f)  # Skip header line
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    gene_id, accession = parts[0], parts[2]
                    if accession not in mapping:
                        mapping[accession] = set()
                    mapping[accession].add(gene_id)
    except FileNotFoundError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    return mapping
def generate_report(json_file_path, output_tsv_path, geneid_accession_path):
    """
    Reads a JSON file from the parse_sam step and generates a final
    tab-separated report with the specified format.
    """

    # Load GeneID to accession mapping if provided
    accession_to_geneid = {}
    if geneid_accession_path:
        accession_to_geneid = load_geneid_accession_map(geneid_accession_path)

    # Read and parse the JSON file
    parsed_data = {}
    try:
        with open(json_file_path, 'r') as f:
            parsed_data = json.load(f)
    except (FileNotFoundError, json.JSONDecodeError) as e:
        print(f"Error reading or parsing JSON file {json_file_path}: {e}", file=sys.stderr)
        sys.exit(1)

    # Open the output TSV file for writing
    with open(output_tsv_path, 'w', newline='') as f_out:
        writer = csv.writer(f_out, delimiter='\t')

        # Write the new, corrected header
        writer.writerow(['#ID', 'mismatch_level', 'num_of_matched_geneids', 'num_of_matched_accessions', 'matched_geneid', 'matched_accession'])

        # Iterate through each oligo in the JSON data
        for oligo_id, data in parsed_data.items():

            # Iterate through each mismatch level for the current oligo
            # The keys are strings in JSON, so we sort them numerically
            sorted_mismatch_levels = sorted(data.get('mismatch_level', {}).keys(), key=int)

            for mismatch_level in sorted_mismatch_levels:
                mismatch_info = data['mismatch_level'][mismatch_level]
                accessions = mismatch_info.get('accessions', [])

                # Map accessions to GeneIDs
                geneids = set()
                if accession_to_geneid:
                    for acc in accessions:
                        if matched_id_set := accession_to_geneid.get(acc):
                            geneids.update(matched_id_set)

                # Prepare the final values for the new columns
                num_of_matched_geneids = len(geneids)
                if num_of_matched_geneids > 10:
                    matched_geneid = 'too_many_to_record'
                else:
                    matched_geneid = ','.join(geneids) if geneids else 'NA'
                    
                num_of_matched_accessions = len(accessions)
                if num_of_matched_accessions > 10:
                    matched_accession = 'too_many_to_record'
                else:
                    matched_accession = ','.join(accessions)

                # Write the final row with the new column order
                writer.writerow([
                    oligo_id,
                    mismatch_level,
                    num_of_matched_geneids,
                    num_of_matched_accessions,
                    matched_geneid,
                    matched_accession
                ])

def main():
    parser = argparse.ArgumentParser(description="Generate a TSV report from a parsed SAM JSON file.")
    parser.add_argument("--json", required=True, help="Input JSON file path.")
    parser.add_argument("--output", required=True, help="Output TSV file path.")
    parser.add_argument("--geneid_accession", required=True, help="Convert accessions to GeneID.")
    args = parser.parse_args()

    generate_report(args.json, args.output, args.geneid_accession)

if __name__ == "__main__":
    main()