process CALCULATE_TARGET_ACCESSIBILITY {
    tag "${params.run_id} - $gene_id - Calculate Target Accessibility"

    input:
    tuple val(gene_id), path(target_gene)

    output:
    path("${gene_id}.target_accessibility.tsv"), emit: target_accessibility

    script:
    def output_accessibility = "${gene_id}.target_accessibility.tsv"

    """
    calculate_target_accessibility.py \
        --input_fasta ${target_gene} \
        --output ${output_accessibility} \
        --winsize ${params.plfold_winsize} \
        --span ${params.plfold_span} \
        --ulength ${params.plfold_ulength} \
        --surrounding_region_length ${params.surrounding_region_length} \
        --oligo_length ${params.oligo_length} \
        --offset_5_prime ${params.offset_5_prime} \
        --gene_id ${gene_id}
    """
}