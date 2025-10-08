#!/usr/bin/env python

import argparse
import sys
import pandas as pd

def load_seqs(seq_file):
    """Loads sequences from a TSV file into a pandas DataFrame."""
    try:
        return pd.read_csv(seq_file, sep='\t')
    except Exception as e:
        print(f"Error loading sequences file: {e}", file=sys.stderr)
        sys.exit(1)

def filter_gc_content(seqs, min_gc, max_gc):
    """Filters sequences based on GC content range."""
    return seqs[(seqs['GC_Content'] >= min_gc) & (seqs['GC_Content'] <= max_gc)]

def filter_microrna_hits(seqs, microrna_hits_threshold):
    """Filters sequences based on microRNA hits threshold."""
    return seqs[seqs['MicroRNA_Hits'] <= microrna_hits_threshold]

def has_forbidden_motif(seq, forbidden_motifs):
    """Checks if a sequence contains any forbidden motifs."""
    return any(motif.upper() in seq for motif in forbidden_motifs)

def filter_forbidden_motifs(seqs, forbidden_motifs):
    """Filters sequences containing any forbidden motifs."""
    return seqs[
        ~seqs['Oligo'].apply(
            lambda x: has_forbidden_motif(x, forbidden_motifs)
        )
    ]

def filter_sequences(seqs, min_gc, max_gc, microrna_hits_threshold, forbidden_motifs):
    """Applies all filters to the seqs DataFrame."""
    seqs = filter_gc_content(seqs, min_gc, max_gc)
    seqs = filter_microrna_hits(seqs, microrna_hits_threshold)
    forbidden_motifs_list = [motif.strip() for motif in forbidden_motifs.split(',') if motif.strip()]
    seqs = filter_forbidden_motifs(seqs, forbidden_motifs_list)
    return seqs

def main():
    parser = argparse.ArgumentParser(description="Filter sequences by GC content, microRNA hits, and forbidden motifs.")
    parser.add_argument("--seq_file", type=str, required=True, help="Path to the sequence TSV file.")
    parser.add_argument("--min_gc", type=float, default=40.0, help="Minimum GC content percentage.")
    parser.add_argument("--max_gc", type=float, default=60.0, help="Maximum GC content percentage.")
    parser.add_argument("--microrna_hits_threshold", type=int, default=1, help="Maximum allowed microRNA hits.")
    parser.add_argument("--forbidden_motifs", type=str, default="", help="Comma-separated list of forbidden motifs.")
    parser.add_argument("--output_file", required=True, help="Output TSV file for filtered sequences.")
    args = parser.parse_args()

    seqs = load_seqs(args.seq_file)
    filtered = filter_sequences(seqs, args.min_gc, args.max_gc, args.microrna_hits_threshold, args.forbidden_motifs)
    filtered.to_csv(args.output_file, sep='\t', index=False)

if __name__ == "__main__":
    main()

