#!/usr/bin/env python

import argparse
import sys
import pandas as pd

def calculate_gc(seq):
    """Calculates the GC content of a DNA sequence."""
    if not seq:
        return 0.0
    gc_count = sum(base in ('G', 'C') for base in seq)
    return (gc_count / len(seq)) * 100

def reverse_complement(seq):
    """Generates the reverse complement of a DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return "".join(complement.get(base, base) for base in reversed(seq.upper()))

def convert_dna_to_rna(seq):
    """Converts a DNA sequence to an RNA sequence by replacing T with U."""
    return seq.upper().replace('T', 'U')

def load_weight_matrix(weight_matrix_file):
    """Loads a weight matrix from a file into a pandas DataFrame."""
    try:
        weight_matrix = pd.read_csv(weight_matrix_file, sep='\t', index_col=0)
    except Exception as e:
        print(f"Error loading weight matrix file: {e}", file=sys.stderr)
        sys.exit(1)
    return weight_matrix

def calc_seq_score(seq, weight_matrix):
    """Calculates a score for the sequence based on a provided weight matrix."""
    seq = convert_dna_to_rna(seq)
    weight_matrix = load_weight_matrix(weight_matrix)
    score = 0.0
    for i, nucleotide in enumerate(seq):
        if nucleotide in weight_matrix.columns:
            score += weight_matrix.loc[i + 1, nucleotide]
        else:
            print(f"Warning: Nucleotide '{nucleotide}' at position {i} not found in weight matrix.", file=sys.stderr)
    return score

def load_microrna_seeds(microrna_seeds_file):
    """Loads microRNA seeds from a file into a set."""
    try:
        with open(microrna_seeds_file, 'r') as f:
            seeds = {line.strip().upper() for line in f if line.strip()}
    except Exception as e:
        print(f"Error loading microRNA seeds file: {e}", file=sys.stderr)
        sys.exit(1)
    return seeds

def calc_microrna_hits(seq, microrna_seeds):
    """Calculates the number of microRNA seed matches in the sequence."""
    seq = convert_dna_to_rna(seq)
    microrna_seeds = load_microrna_seeds(microrna_seeds)

    hit_count = 0
    seed_length = len(next(iter(microrna_seeds))) if microrna_seeds else 0

    for i in range(len(seq) - seed_length + 1):
        subseq = seq[i:i + seed_length]
        if subseq in microrna_seeds:
            hit_count += 1

    return hit_count

def load_cds_regions(cds_region_file):
    """Loads CDS regions from a file into a dictionary."""
    cds_regions = {}
    try:
        with open(cds_region_file, 'r') as f:
            next(f)  # Skip header
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) == 3:
                    accession, start, end = parts
                    cds_regions[accession] = (int(start), int(end))
    except Exception as e:
        print(f"Error loading CDS regions file: {e}", file=sys.stderr)
        sys.exit(1)
    return cds_regions

def generate_sequences(input_fasta, output, surrounding_region_length, 
                    offset_5_prime, oligo_length, offset_refseq_seed, refseq_seed_length, 
                    offset_microrna, microrna_seed_length, weight_matrix, microrna_seeds, cds_region_file):
    # sourcery skip: low-code-quality
    """
    Reads a FASTA file, extracts the surrounding region of specified length,
    and writes it to an output FASTA file.
    """
    # Load CDS regions
    cds_regions = load_cds_regions(cds_region_file)
    
    # --- Read the input FASTA file ---
    sequence = ""
    accession = ""
    with open(input_fasta, 'r') as f_in:
        for line in f_in:
            if line.startswith('>'):
                accession = line[1:].strip().split()[0]
                if accession not in cds_regions:
                    print(f"Error: Accession {accession} not found in CDS regions file.", file=sys.stderr)
                    sys.exit(1)
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
    header = "#ID\tSurrounding_Region\tOligo\tRegion\tGC_Content\tRefseq_Seed\tOligo_RC\tMicroRNA_Seed\tMicroRNA_Hits\tScore\n"
    with open(output, 'w') as f_out:
        f_out.write(header)
        for i in range(end):
            # Extract the surrounding region
            surrounding_region = sequence[i:i + surrounding_region_length].upper()
            score = calc_seq_score(surrounding_region, weight_matrix)
            
            # Extract the oligo and other oligo based on sequence
            oligo = surrounding_region[offset_5_prime:offset_5_prime + oligo_length]
            gc_oligo = calculate_gc(oligo)
            refseq_seed = oligo[offset_refseq_seed:offset_refseq_seed + refseq_seed_length]
            
            # Identify if the oligo overlaps with the CDS region
            region = ""
            cds_start, cds_end = cds_regions[accession]
            oligo_start = i + offset_5_prime + 1  # 1-based
            oligo_end = oligo_start + oligo_length - 1
            # Determine the region type
            if oligo_end < cds_start: # completely before CDS
                region = "5UTR"
            elif oligo_start < cds_start: # overlaps 5' UTR and CDS
                region = "5UTR_CDS"
            elif oligo_start <= cds_end and oligo_end > cds_end: # overlaps CDS and 3' UTR
                region = "CDS_3UTR"
            elif oligo_start > cds_end: # completely after CDS
                region = "3UTR"
            else: # completely within CDS
                region = "CDS"
            
            # Generate the reverse complement of the oligo and microRNA seed
            oligo_rc = reverse_complement(oligo)
            microrna_seed = convert_dna_to_rna(oligo_rc[offset_microrna:offset_microrna + microrna_seed_length])
            microrna_hits = calc_microrna_hits(microrna_seed, microrna_seeds)
            
            result_data = [
                i,
                surrounding_region,
                oligo,
                region,
                f"{gc_oligo:.2f}",
                refseq_seed,
                oligo_rc,
                microrna_seed,
                microrna_hits,
                score
            ]

            output_str = "\t".join(map(str, result_data)) + "\n"
            # Write the output string to the file
            f_out.write(output_str)

def main():
    parser = argparse.ArgumentParser(description="Generate seqsuences and metadata from a FASTA file.")
    parser.add_argument("--input_fasta", required=True, help="Input FASTA file")
    parser.add_argument("--output", required=True, help="Output metadata file")
    parser.add_argument("--surrounding_region_length", type=int, required=True, help="Length of surrounding region")
    parser.add_argument("--offset_5_prime", type=int, required=True, help="5' offset")
    parser.add_argument("--oligo_length", type=int, required=True, help="Oligo length")
    parser.add_argument("--offset_refseq_seed", type=int, required=True, help="Refseq seed offset")
    parser.add_argument("--refseq_seed_length", type=int, required=True, help="Refseq seed length")
    parser.add_argument("--offset_microrna", type=int, required=True, help="MicroRNA offset")
    parser.add_argument("--microrna_seed_length", type=int, required=True, help="MicroRNA seed length")
    parser.add_argument("--weight_matrix", type=str, required=True, help="Weight matrix file")
    parser.add_argument("--microrna_seeds", type=str, required=True, help="MicroRNA seeds file")
    parser.add_argument("--cds_region", type=str, required=True, help="CDS region file")
    args = parser.parse_args()

    generate_sequences(
        args.input_fasta, 
        args.output, 
        args.surrounding_region_length,
        args.offset_5_prime,
        args.oligo_length,
        args.offset_refseq_seed,
        args.refseq_seed_length,
        args.offset_microrna,
        args.microrna_seed_length,
        args.weight_matrix,
        args.microrna_seeds,
        args.cds_region
    )

if __name__ == "__main__":
    main()
    

