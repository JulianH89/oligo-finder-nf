#!/usr/bin/env python

import argparse
import sys
import textwrap

def generate_oligos(input_fasta, output_fasta, oligo_length, offset_5, offset_3):
    """
    Reads a single-sequence FASTA file, applies 5' and 3' offsets,
    and generates oligo candidates of a specified length.
    """
    # --- Read the input FASTA file ---
    sequence = ""
    with open(input_fasta, 'r') as f_in:
        for line in f_in:
            if line.startswith('>'):
                continue
            sequence += line.strip()

    if not sequence:
        print(f"Error: No sequence found in {input_fasta}", file=sys.stderr)
        sys.exit(1)

    # --- Apply offsets to trim the sequence ---
    if offset_3 == 0:
        trimmed_sequence = sequence[offset_5:]
    else:
        trimmed_sequence = sequence[offset_5:-offset_3]

    if len(trimmed_sequence) < oligo_length:
        print(f"Error: Trimmed sequence length ({len(trimmed_sequence)}) is less than the oligo length ({oligo_length}).", file=sys.stderr)
        sys.exit(1)

    # --- Generate oligo candidates and write to output FASTA ---
    with open(output_fasta, 'w') as f_out:
        for i in range(len(trimmed_sequence) - oligo_length + 1):
            oligo_seq = trimmed_sequence[i : i + oligo_length]
            oligo_id = f"oligo_{i + 1}"
            f_out.write(f">{oligo_id}\n")
            f_out.write("\n".join(textwrap.wrap(oligo_seq, 60)))
            f_out.write("\n")

def main():
    parser = argparse.ArgumentParser(description="Slice a gene sequence into oligo candidates.")
    parser.add_argument("--input-fasta", required=True, help="Input FASTA file containing the target gene.")
    parser.add_argument("--output-fasta", required=True, help="Output FASTA file for oligo candidates.")
    parser.add_argument("--oligo-length", required=True, type=int, help="The length of each oligo.")
    parser.add_argument("--offset-5", required=True, type=int, help="Number of bases to trim from the 5' end.")
    parser.add_argument("--offset-3", required=True, type=int, help="Number of bases to trim from the 3' end.")
    args = parser.parse_args()

    generate_oligos(args.input_fasta, args.output_fasta, args.oligo_length, args.offset_5, args.offset_3)

if __name__ == "__main__":
    main()