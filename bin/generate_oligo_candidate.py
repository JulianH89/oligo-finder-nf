#!/usr/bin/env python

import argparse
import sys
from collections import OrderedDict

def calculate_gc(seq):
    """Calculates the GC content of a DNA sequence."""
    if not seq:
        return 0.0
    gc_count = seq.upper().count('G') + seq.upper().count('C')
    return (gc_count / len(seq)) * 100

def reverse_complement(seq):
    """Generates the reverse complement of a DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return "".join(complement.get(base, base) for base in reversed(seq.upper()))

def has_forbidden_motif(seq, motifs):
    """Checks if a sequence contains any forbidden motifs."""
    upper_seq = seq.upper()
    return any(motif.upper() in upper_seq for motif in motifs)

def generate_and_filter_oligos(input_fasta, output_fasta, oligo_length, offset_5, offset_3, min_gc, max_gc, forbidden_motifs):
    """
    Reads a FASTA file, generates oligo candidates, filters them,
    and writes the valid candidates to an output FASTA file.
    """
    # --- Read the input FASTA file ---
    sequence = ""
    with open(input_fasta, 'r') as f_in:
        for line in f_in:
            if not line.startswith('>'):
                sequence += line.strip()

    if not sequence:
        print(f"Error: No sequence found in {input_fasta}", file=sys.stderr)
        sys.exit(1)

    # --- Apply offsets ---
    if offset_3 == 0:
        trimmed_sequence = sequence[offset_5:]
    else:
        trimmed_sequence = sequence[offset_5:-offset_3]

    if len(trimmed_sequence) < oligo_length:
        print(f"Error: Trimmed sequence length ({len(trimmed_sequence)}) is less than the oligo length ({oligo_length}).", file=sys.stderr)
        sys.exit(1)

    # --- Generate, filter, and write oligo candidates ---
    # This list comprehension now filters out any empty strings from the motifs list.
    forbidden_motifs_list = [motif.strip() for motif in forbidden_motifs.split(',') if motif.strip()]
    passed_candidates = OrderedDict()
    oligo_count = 0

    end = len(trimmed_sequence) - oligo_length + 1
    for i in range(end):
        fwd_seq = trimmed_sequence[i : i + oligo_length].upper()

        # --- Filter 1: Uniqueness ---
        if fwd_seq in passed_candidates:
            continue

        rev_comp_seq = reverse_complement(fwd_seq)

        # --- Filter 2: GC Content ---
        gc_fwd = calculate_gc(fwd_seq)
        # if not (min_gc <= gc_fwd <= max_gc):
        #     continue

        gc_rev = calculate_gc(rev_comp_seq)
        if not (min_gc <= gc_rev <= max_gc):
            continue

        # --- Filter 3: Forbidden Motifs ---
        # if has_forbidden_motif(fwd_seq, forbidden_motifs_list):
        #     continue
        if has_forbidden_motif(rev_comp_seq, forbidden_motifs_list):
            continue

        # If all checks pass, it's a valid candidate
        passed_candidates[fwd_seq] = None
        oligo_count += 1
    
    # --- Write final valid candidates to file ---
    with open(output_fasta, 'w') as f_out:
        for i, seq in enumerate(passed_candidates.keys(), 1):
            f_out.write(f">oligo_{i}\n")
            f_out.write(f"{seq}\n")

def main():
    parser = argparse.ArgumentParser(description="Generate and filter oligo candidates from a gene sequence.")
    parser.add_argument("--input-fasta", required=True, help="Input FASTA file containing the target gene.")
    parser.add_argument("--output-fasta", required=True, help="Output FASTA file for valid oligo candidates.")
    parser.add_argument("--oligo-length", required=True, type=int, help="The length of each oligo.")
    parser.add_argument("--offset-5", required=True, type=int, help="Number of bases to trim from the 5' end.")
    parser.add_argument("--offset-3", required=True, type=int, help="Number of bases to trim from the 3' end.")
    parser.add_argument("--min-gc", required=True, type=float, help="Minimum GC content percentage.")
    parser.add_argument("--max-gc", required=True, type=float, help="Maximum GC content percentage.")
    parser.add_argument("--forbidden-motifs", required=True, help="Comma-separated list of forbidden motifs.")
    args = parser.parse_args()

    generate_and_filter_oligos(
        args.input_fasta, 
        args.output_fasta, 
        args.oligo_length, 
        args.offset_5, 
        args.offset_3, 
        args.min_gc, 
        args.max_gc, 
        args.forbidden_motifs
    )

if __name__ == "__main__":
    main()