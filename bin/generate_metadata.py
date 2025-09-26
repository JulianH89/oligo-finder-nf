#!/usr/bin/env python

import argparse
import sys
import pandas as pd

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

def convert_dna_to_rna(seq):
    """Converts a DNA sequence to an RNA sequence by replacing T with U."""
    return seq.upper().replace('T', 'U')

def calc_seq_score(seq, weight_matrix):
    """Calculates a score for the sequence based on a provided weight matrix."""
    seq = convert_dna_to_rna(seq)
    weight_matrix = pd.read_csv(weight_matrix, sep = '\t', index_col=0)
    score = 0.0
    for i, nucleotide in enumerate(seq):
        if nucleotide in weight_matrix.columns:
            score += weight_matrix.loc[i + 1, nucleotide]
        else:
            print(f"Warning: Nucleotide '{nucleotide}' at position {i} not found in weight matrix.", file=sys.stderr)
    return score

def generate_surrounding_region(input_fasta, output, surrounding_region_length, 
                    offset_5_prime, oligo_length, offset_gene_region, gene_region_length, 
                    offset_microrna, microrna_seed_length, weight_matrix):
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

    # --- Validate lengths ---
    if len(sequence) < surrounding_region_length:
        print(f"Error: Sequence length ({len(sequence)}) is less than the surrounding region length ({surrounding_region_length}).", file=sys.stderr)
        sys.exit(1)

    end = len(sequence) - surrounding_region_length + 1
    
    # --- Write to output metadata file ---
    header = "#ID\tSurrounding_Region\tOligo\tGC_Content\tGene_Region\tOligo_RC\tMicroRNA_Seed\tScore\n"
    with open(output, 'w') as f_out:
        f_out.write(header)
        for i in range(end):
            # Extract the surrounding region
            surrounding_region = sequence[i:i + surrounding_region_length].upper()
            score = calc_seq_score(surrounding_region, weight_matrix)
            
            # Extract the oligo and other oligo based on sequence
            oligo = surrounding_region[offset_5_prime:offset_5_prime + oligo_length]
            gc_oligo = calculate_gc(oligo)
            gene_region = oligo[offset_gene_region:offset_gene_region + gene_region_length]
            
            # Generate the reverse complement of the oligo and microRNA seed
            oligo_rc = reverse_complement(oligo)
            microrna_seed = convert_dna_to_rna(oligo_rc[offset_microrna:offset_microrna + microrna_seed_length])
            
            result_data = [
                i,
                surrounding_region,
                oligo,
                f"{gc_oligo:.2f}",
                gene_region,
                oligo_rc,
                microrna_seed,
                score
            ]

            output_str = "\t".join(map(str, result_data)) + "\n"
            # Write the output string to the file
            f_out.write(output_str)

def main():
    parser = argparse.ArgumentParser(description="Generate surrounding region from FASTA")
    parser.add_argument("--input-fasta", required=True, help="Input FASTA file")
    parser.add_argument("--output", required=True, help="Output metadata file")
    parser.add_argument("--surrounding-region-length", type=int, required=True, help="Length of surrounding region")
    parser.add_argument("--offset_5_prime", type=int, required=True, help="5' offset")
    parser.add_argument("--oligo-length", type=int, required=True, help="Oligo length")
    parser.add_argument("--offset_gene_region", type=int, required=True, help="Gene region offset")
    parser.add_argument("--gene_region_length", type=int, required=True, help="Gene region length")
    parser.add_argument("--offset_microrna", type=int, required=True, help="MicroRNA offset")
    parser.add_argument("--microrna_seed_length", type=int, required=True, help="MicroRNA seed length")
    parser.add_argument("--weight_matrix", type=str, help="Weight matrix file (optional)")
    args = parser.parse_args()

    generate_surrounding_region(
        args.input_fasta, 
        args.output, 
        args.surrounding_region_length,
        args.offset_5_prime,
        args.oligo_length,
        args.offset_gene_region,
        args.gene_region_length,
        args.offset_microrna,
        args.microrna_seed_length,
        args.weight_matrix
    )

if __name__ == "__main__":
    main()
    

