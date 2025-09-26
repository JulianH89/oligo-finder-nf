process GENERATE_SURROUNDING_REGION {
    tag "${params.run_id} - $gene_id - Generate Surrounding Region"

    container 'python:3.10'

    input:
    tuple val(gene_id), path(target_gene)

    output:
    tuple val(gene_id), path("${target_gene.baseName}.surrounding.fa"), emit: surrounding_fasta

    script:
    def output_fasta = "${target_gene.baseName}.surrounding.fa"
    """
    generate_metadata.py \\
        --input-fasta ${target_gene} \\
        --output-fasta ${output_fasta} \\
        --surrounding-region-length ${params.surrounding_region_length} \\
        --offset_5_prime ${params.offset_5_prime}
        --oligo-length ${params.oligo_length}
    """
}