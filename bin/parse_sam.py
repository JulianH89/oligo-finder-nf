#!/usr/bin/env python

import argparse
import json
import sys
from collections import defaultdict

def calculate_gc_content(sequence):
    """Calculates the GC content of a DNA sequence."""
    if not sequence:
        return 0.0
    gc_count = sequence.upper().count('G') + sequence.upper().count('C')
    return (gc_count / len(sequence)) * 100

def parse_sam_file(sam_file_path):
    """
    Parses a Bowtie SAM file and structures the data as a nested dictionary.
    """
    oligos = {}

    with open(sam_file_path, 'r') as f:
        for line in f:
            if line.startswith('@'):
                continue

            fields = line.strip().split('\t')
            
            if len(fields) < 11:
                continue

            oligo_id = fields[0]
            accession = fields[2]
            sequence = fields[9]
            
            # --- OPTIMIZATION 1: Faster Tag Lookup ---
            # Instead of looping, create a dictionary of the optional tags for instant lookup.
            tags = {tag.split(':')[0]: tag.split(':')[-1] for tag in fields[11:]}
            
            num_mismatches_str = tags.get('NM')
            if not num_mismatches_str:
                continue
            
            try:
                num_mismatches = int(num_mismatches_str)
            except ValueError:
                continue

            if oligo_id not in oligos:
                # The redundant 'oligo_id' field has been removed from this dictionary
                oligos[oligo_id] = {
                    'sequence': sequence,
                    'gc_content': round(calculate_gc_content(sequence), 2),
                    # --- OPTIMIZATION 2 (MAJOR): Use a set instead of a list ---
                    # This makes checking for existing accessions much faster.
                    'mismatches': defaultdict(lambda: {'accessions': set()})
                }

            # Adding to a set is extremely fast, even for thousands of items.
            oligos[oligo_id]['mismatches'][num_mismatches]['accessions'].add(accession)

    # Convert sets back to lists for JSON compatibility and convert defaultdicts
    for oligo_id, data in oligos.items():
        converted_mismatches = {}
        for mismatch_level, mismatch_data in data['mismatches'].items():
            converted_mismatches[mismatch_level] = {
                'accessions': list(mismatch_data['accessions'])
            }
        data['mismatches'] = converted_mismatches

    return oligos

def main():
    parser = argparse.ArgumentParser(description="Parse a SAM file to a structured JSON format.")
    parser.add_argument("--sam", required=True, help="Input SAM file path.")
    parser.add_argument("--output", required=True, help="Output JSON file path.")
    args = parser.parse_args()

    try:
        parsed_data = parse_sam_file(args.sam)
        with open(args.output, 'w') as f_out:
            json.dump(parsed_data, f_out, indent=4)
    except Exception as e:
        print(f"Error processing file {args.sam}: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()