process GENERATE_METADATA {
    tag "${params.run_id} - $gene_id - Generate Metadata"
    publishDir "${params.outdir}/${params.run_id}/${gene_id}", mode: 'copy'

    input:
    tuple val(gene_id), path(target_gene)

    output:
    tuple val(gene_id), path("${gene_id}.metadata.tsv"), emit: seq_metadata

    script:
    def seq_metadata = "${gene_id}.metadata.tsv"
    """
    generate_metadata.py \\
        --input_fasta ${target_gene} \\
        --output ${seq_metadata} \\
        --surrounding_region_length ${params.surrounding_region_length} \\
        --offset_5_prime ${params.offset_5_prime} \\
        --oligo_length ${params.oligo_length} \\
        --offset_refseq_seed ${params.offset_refseq_seed} \\
        --refseq_seed_length ${params.refseq_seed_length} \\
        --offset_microrna ${params.offset_microrna} \\
        --microrna_seed_length ${params.microrna_seed_length} \\
        --weight_matrix ${params.weight_matrix} \\
        --microrna_seeds ${params.microrna_seeds}
    """
}