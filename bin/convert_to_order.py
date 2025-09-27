#!/usr/bin/env python
import sys
import argparse
import pandas as pd

def order_oligo_sense_no_tripurine(oligo, sense_length):
    """
    Converts a sense strand to TriLink format with a fixed modification pattern,
    without the special tripurine handling.

    Args:
        oligo: The input oligonucleotide sequence.
        sense_length: The length of the sense strand to process.

    Returns:
        The formatted string for synthesis.
    """
    sense = oligo[5 : 5 + sense_length]
    mod_nuc_parts = []
    
    # Define modifications for different positions
    full_m_mod = {'C': 'mC', 'U': 'mU', 'A': 'mA', 'G': 'mG'}
    partial_m_mod = {'C': 'mC', 'U': 'mU', 'A': 'A', 'G': 'G'}

    for i, base in enumerate(sense):
        if i < 2:
            mod_nuc_parts.append(full_m_mod.get(base, base) + ".")
        elif 2 <= i < 12:
            mod_nuc_parts.append(partial_m_mod.get(base, base) + ".")
        elif i == 12:
            mod_nuc_parts.append(partial_m_mod.get(base, base) + "#")
        else: # i > 12
            mod_nuc_parts.append(full_m_mod.get(base, base))
            
    return "".join(mod_nuc_parts) + "#mA1"


def order_oligo_sense(oligo, sense_length):
    """
    Converts a sense strand to TriLink format, applying 2'-O-methyl
    modifications to the 3rd consecutive purine (A/G) to avoid
    an immune response.

    Args:
        oligo: The input oligonucleotide sequence.
        sense_length: The length of the sense strand to process.

    Returns:
        The formatted string for synthesis.
    """
    sense = oligo[5 : 5 + sense_length]
    mod_nuc_parts = []
    purine_flag = 0
    
    full_m_mod = {'C': 'mC', 'U': 'mU', 'A': 'mA', 'G': 'mG'}

    for i, base in enumerate(sense):
        mod_base = base
        if i < 2:
            mod_base = full_m_mod.get(base, base)
            mod_nuc_parts.append(mod_base + ".")
        elif 2 <= i <= 12:
            # Tripurine logic
            if base in ('A', 'G'):
                purine_flag += 1
            else:
                purine_flag = 0
            
            purine_modulus = purine_flag % 3
            
            # Apply modifications
            if base in ('C', 'U'):
                mod_base = 'm' + base
            elif base in ('A', 'G') and purine_flag > 0 and purine_modulus == 0:
                mod_base = 'm' + base
            
            linkage = "." if i < 12 else "#"
            mod_nuc_parts.append(mod_base + linkage)
        else: # i > 12
            mod_base = full_m_mod.get(base, base)
            mod_nuc_parts.append(mod_base)

    return "".join(mod_nuc_parts) + "#mA1"


def order_oligo_antisense(revcomp, antisense_length):
    """
    Converts an antisense strand to a synthesis format, applying 2'-Fluoro
    mods and the tripurine (3rd consecutive A/G) 2'-O-methyl mod.

    Args:
        revcomp: The reverse-complemented oligonucleotide sequence.
        antisense_length: The length of the antisense strand to process.

    Returns:
        The formatted string for synthesis.
    """
    antisense = revcomp[1 : 1 + antisense_length]
    mod_nuc_parts = []
    purine_flag = 0
    
    for i, base in enumerate(antisense):
        mod_base = base
        # Tripurine logic applies to the whole sequence
        if base in ('A', 'G'):
            purine_flag += 1
        else:
            purine_flag = 0
        
        purine_modulus = purine_flag % 3
        
        # Apply base modifications
        if base in ('C', 'U'):
            mod_base = 'f' + base
        elif base in ('A', 'G') and purine_flag > 0 and purine_modulus == 0:
            mod_base = 'm' + base
        
        # Determine linkage
        if i < 12:
            mod_nuc_parts.append(mod_base + ".")
        elif 12 <= i < 18:
            mod_nuc_parts.append(mod_base + "#")
        else: # i >= 18
            mod_nuc_parts.append(base) # Original Perl code shows no mods here

    return "PmU." + "".join(mod_nuc_parts)


def order_oligo_sense_fm(oligo, sense_length):
    """
    Converts a sense strand using an alternating 2'-Fluoro / 2'-O-methyl
    modification pattern (FM format).

    Args:
        oligo: The input oligonucleotide sequence.
        sense_length: The length of the sense strand to process.

    Returns:
        The formatted string for synthesis.
    """
    sense = oligo[5 : 5 + sense_length]
    mod_nuc_parts = []
    
    for i, base in enumerate(sense):
        # Determine modification based on position (even/odd)
        mod_prefix = 'f' if i % 2 == 0 else 'm'
        mod_base = mod_prefix + base
        
        # Determine linkage
        if i < (sense_length - 2):
            mod_nuc_parts.append(mod_base + ".")
        elif i == (sense_length - 2):
            mod_nuc_parts.append(mod_base)
        else: # i == (sense_length - 1)
            mod_nuc_parts.append("#" + mod_base)
            
    return "".join(mod_nuc_parts) + "#fA1"


def order_oligo_antisense_fm(revcomp, antisense_length):
    """
    Converts an antisense strand using an alternating 2'-Fluoro / 2'-O-methyl
    modification pattern (FM format).

    Args:
        revcomp: The reverse-complemented oligonucleotide sequence.
        antisense_length: The length of the antisense strand to process.

    Returns:
        The formatted string for synthesis.
    """
    antisense = revcomp[1 : 1 + antisense_length]
    mod_nuc_parts = []
    
    for i, base in enumerate(antisense):
        mod_prefix = 'f' if i % 2 == 0 else 'm'
        mod_base = mod_prefix + base

        # Determine linkage
        if i < (antisense_length - 7):
            mod_nuc_parts.append(mod_base + ".")
        elif i == (antisense_length - 7):
            mod_nuc_parts.append(mod_base)
        else: # i > (antisense_length - 7)
            mod_nuc_parts.append("#" + mod_base)
            
    return "PmU." + "".join(mod_nuc_parts)

def convert_to_order(report_tsv, sense_length, antisense_length, output_tsv):
    # Load the report TSV
    try:
        df = pd.read_csv(report_tsv, sep="\t")
    except Exception as e:
        print(f"Error loading file {report_tsv}: {e}", file=sys.stderr)
        sys.exit(1)

    # Apply the conversion functions to each row
    df['Sense_Tripurine'] = df['Oligo'].replace("U", "T").apply(lambda x: order_oligo_sense(x, sense_length))
    df['Antisense_Tripurine'] = df['Oligo_RC'].replace("U", "T").apply(lambda x: order_oligo_antisense(x, antisense_length))
    df['Sense_FM'] = df['Oligo'].replace("U", "T").apply(lambda x: order_oligo_sense_fm(x, sense_length))
    df['Antisense_FM'] = df['Oligo_RC'].replace("U", "T").apply(lambda x: order_oligo_antisense_fm(x, antisense_length))
    
    new_column_order = [
    '#ID',
    'Surrounding_Region',
    'Oligo',
    'Oligo_RC',
    'GC_Content',
    'Sense_Tripurine',
    'Antisense_Tripurine',
    'Sense_FM',
    'Antisense_FM',
    'Refseq_Seed',
    'Score',
    'mismatch_level',
    'num_of_matched_accessions',
    'MicroRNA_Seed',
    'MicroRNA_Hits',
    'matched_accession'
    ]
    df = df[new_column_order]

    # Save the updated DataFrame to the output TSV
    try:
        df.to_csv(output_tsv, sep="\t", index=False)
    except Exception as e:
        print(f"Error saving file {output_tsv}: {e}", file=sys.stderr)
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description="Convert oligo sequences to TriLink order format.")
    parser.add_argument("--report_tsv", required=True, help="Path to the input report TSV file.")
    parser.add_argument("--sense_length", type=int, default=14, help="Length of the sense strand.")
    parser.add_argument("--antisense_length", type=int, default=19, help="Length of the antisense strand.")
    parser.add_argument("--output_tsv", required=True, help="Path to the output TSV file.")
    args = parser.parse_args()
    
    convert_to_order(args.report_tsv, args.sense_length, args.antisense_length, args.output_tsv)


if __name__ == '__main__':
    main()