#!/usr/bin/env python

import argparse
import sys

def calculate_gc(seq):
    """Calculates the GC content of a DNA sequence."""
    if not seq:
        return 0.0
    gc_count = seq.upper().count('G') + seq.upper().count('C')
    return (gc_count / len(seq)) * 100

def generate_surrounding_region(input_fasta, output_fasta, surrounding_region_length, 
                    offset_5_prime, oligo_length):
    """
    Reads a FASTA file, extracts the surrounding region of specified length,
    and writes it to an output FASTA file.
    """
    # --- Read the input FASTA file ---
    sequence = ""
    header = ""
    with open(input_fasta, 'r') as f_in:
        for line in f_in:
            if line.startswith('>'):
                header = line.strip()
            else:
                sequence += line.strip()

    if not sequence:
        print(f"Error: No sequence found in {input_fasta}", file=sys.stderr)
        sys.exit(1)

    # --- Extract surrounding region ---
    if len(sequence) < surrounding_region_length:
        print(f"Error: Sequence length ({len(sequence)}) is less than the surrounding region length ({surrounding_region_length}).", file=sys.stderr)
        sys.exit(1)

    end = len(sequence) - surrounding_region_length + 1
    # --- Write to output FASTA file ---
    with open(output_fasta, 'w') as f_out:
        for i in range(end):
            surrounding_region = sequence[i:i + surrounding_region_length].upper()
            oligo = surrounding_region[offset_5_prime:offset_5_prime + oligo_length]
            gc_oligo = calculate_gc(oligo)

            output_str = f"{i}\t{surrounding_region}\t{oligo}\t{gc_oligo:<.2f}\n"
            # Write the output string to the file
            f_out.write(output_str)

def main():
    parser = argparse.ArgumentParser(description="Generate surrounding region from FASTA")
    parser.add_argument("--input-fasta", required=True, help="Input FASTA file")
    parser.add_argument("--output-fasta", required=True, help="Output FASTA file")
    parser.add_argument("--surrounding-region-length", type=int, required=True, help="Length of surrounding region")
    parser.add_argument("--offset_5_prime", type=int, required=True, help="5' offset")
    parser.add_argument("--oligo-length", type=int, required=True, help="Oligo length")
    args = parser.parse_args()

    generate_surrounding_region(
        args.input_fasta, 
        args.output_fasta, 
        args.surrounding_region_length,
        args.offset_5_prime,
        args.oligo_length
        )

if __name__ == "__main__":
    main()
    

