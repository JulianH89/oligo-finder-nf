#!/usr/bin/env python

import argparse
import sys
import pandas as pd

def load_metadata(metadata_file):
    """Loads metadata from a TSV file into a pandas DataFrame."""
    try:
        metadata = pd.read_csv(metadata_file, sep='\t')
        return metadata
    except Exception as e:
        print(f"Error loading metadata file: {e}", file=sys.stderr)
        sys.exit(1)

def filter_gc_content(metadata, min_gc, max_gc):
    """Filters sequences based on GC content range."""
    filtered = metadata[(metadata['GC_Content'] >= min_gc) & (metadata['GC_Content'] <= max_gc)]
    return filtered

def filter_microrna_hits(metadata, microrna_hits_threshold):
    """Filters sequences based on microRNA hits threshold."""
    filtered = metadata[metadata['MicroRNA_Hits'] <= microrna_hits_threshold]
    return filtered
 
def has_forbidden_motif(seq, forbidden_motifs):
    """Checks if a sequence contains any forbidden motifs."""
    return any(motif.upper() in seq for motif in forbidden_motifs)

def filter_forbidden_motifs(metadata, forbidden_motifs):
    """Filters sequences containing any forbidden motifs."""
    filtered = metadata[~metadata['Oligo'].apply(lambda x: has_forbidden_motif(x, forbidden_motifs))]
    return filtered

def filter_sequences(metadata, min_gc, max_gc, microrna_hits_threshold, forbidden_motifs):
    """Applies all filters to the metadata DataFrame."""
    metadata = filter_gc_content(metadata, min_gc, max_gc)
    metadata = filter_microrna_hits(metadata, microrna_hits_threshold)
    forbidden_motifs_list = [motif.strip() for motif in forbidden_motifs.split(',') if motif.strip()]
    metadata = filter_forbidden_motifs(metadata, forbidden_motifs_list)
    return metadata

def main():
    parser = argparse.ArgumentParser(description="Filter oligo candidates by GC content, microRNA hits, and forbidden motifs.")
    parser.add_argument("--metadata_file", type=str, required=True, help="Path to the metadata TSV file.")
    parser.add_argument("--min_gc", type=float, default=40.0, help="Minimum GC content percentage.")
    parser.add_argument("--max_gc", type=float, default=60.0, help="Maximum GC content percentage.")
    parser.add_argument("--microrna_hits_threshold", type=int, default=1, help="Maximum allowed microRNA hits.")
    parser.add_argument("--forbidden_motifs", type=str, default="", help="Comma-separated list of forbidden motifs.")
    parser.add_argument("--output_file", required=True, help="Output TSV file for filtered oligos.")
    args = parser.parse_args()

    metadata = load_metadata(args.metadata_file)
    filtered = filter_sequences(metadata, args.min_gc, args.max_gc, args.microrna_hits_threshold, args.forbidden_motifs)
    filtered.to_csv(args.output_file, sep='\t', index=False)

if __name__ == "__main__":
    main()

