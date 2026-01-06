#!/usr/bin/env python

import RNA
import argparse

def load_sequence(input_fasta):
    """Loads a RNA sequence."""
    with open(input_fasta, "r") as f:
        lines = f.readlines()
        seq = "".join(
            line.strip().replace("T", "U")
            for line in lines
            if not line.startswith(">")
        )
    return seq.upper()


def calculate_accessibility(gene_id, input_fasta, output, winsize, span, ulength, surrounding_region_length, oligo_length, offset_5_prime):
    # sourcery skip: avoid-builtin-shadow
    """Calculates the target accessibility of a RNA sequence using RNAplfold."""
    
    ## Load sequence
    seq = load_sequence(input_fasta)
    
    ## Run RNAplfold
    pl_matrix = RNA.pfl_fold_up(seq, ulength, winsize, span)
    
    ## Parse and print accessibility results
    results = []
    for i in range(ulength, len(pl_matrix)):
        if i >= ulength:
            accessibility = pl_matrix[i][ulength]
            results.append(accessibility)
            
    results = results[offset_5_prime:-(surrounding_region_length - oligo_length - offset_5_prime)]
    
    with open(output, "w") as out_f:
        out_f.write("#ID\tAccessibility\n")
        for i in range(len(results)):
            id = f"{gene_id}_{i+1}"
            accessibility = results[i]
            out_f.write(f"{id}\t{accessibility:.6f}\n")
    


def main():
    parser = argparse.ArgumentParser(description="Calculate target accessibility of RNA sequences.")
    parser.add_argument("--gene_id", type=str, help="Gene ID/accession.")
    parser.add_argument("--input_fasta", type=str, help="target gene sequence to analyze.")
    parser.add_argument("--output", type=str, help="Output file to save accessibility results.")
    parser.add_argument("--winsize", type=int, default=70, help="Window size for RNAplfold.")
    parser.add_argument("--span", type=int, default=50, help="Span for RNAplfold.")
    parser.add_argument("--ulength", type=int, default=20, help="Maximum length of unpaired regions.")
    parser.add_argument("--surrounding_region_length", type=int, help="Length of surrounding region.")
    parser.add_argument("--oligo_length", type=int, help="Length of oligo.")
    parser.add_argument("--offset_5_prime", type=int, help="Offset for 5' end of oligo.")

    args = parser.parse_args()

    calculate_accessibility(args.gene_id, args.input_fasta, args.output, args.winsize, args.span, args.ulength, args.surrounding_region_length, args.oligo_length, args.offset_5_prime)

if __name__ == "__main__":
    main()