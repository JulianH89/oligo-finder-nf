#!/usr/bin/env python

import RNA
import argparse
import sys


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


def calculate_accessibility(input_fasta, output, winsize, span, ulength):
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
            results.append((i - ulength + 1, accessibility))

    with open(output, "w") as out_f:
        for position, accessibility in results:
            out_f.write(f"{position}\t{accessibility:.6f}\n")
    


def main():
    parser = argparse.ArgumentParser(description="Calculate target accessibility of RNA sequences.")
    parser.add_argument("--input_fasta", type=str, help="target gene sequence to analyze.")
    parser.add_argument("--output", type=str, help="Output file to save accessibility results.")
    parser.add_argument("--winsize", type=int, default=70, help="Window size for RNAplfold.")
    parser.add_argument("--span", type=int, default=50, help="Span for RNAplfold.")
    parser.add_argument("--ulength", type=int, default=20, help="Maximum length of unpaired regions.")

    args = parser.parse_args()

    calculate_accessibility(args.input_fasta, args.output, args.winsize, args.span, args.ulength)

if __name__ == "__main__":
    main()