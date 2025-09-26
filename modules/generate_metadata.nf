process GENERATE_METADATA {
    tag "${params.run_id} - $gene_id - Generate Metadata"
    publishDir "${params.outdir}/${params.run_id}/${gene_id}", mode: 'copy'

    input:
    tuple val(gene_id), path(target_gene)

    output:
    tuple val(gene_id), path("${gene_id}.metadata.txt"), emit: seq_metadata

    script:
    def seq_metadata = "${gene_id}.metadata.txt"
    """
    generate_metadata.py \\
        --input-fasta ${target_gene} \\
        --output ${seq_metadata} \\
        --surrounding-region-length ${params.surrounding_region_length} \\
        --offset_5_prime ${params.offset_5_prime} \\
        --oligo-length ${params.oligo_length} \\
        --offset_gene_region ${params.offset_gene_region} \\
        --gene_region_length ${params.gene_region_length} \\
        --offset_microrna ${params.offset_microrna} \\
        --microrna_seed_length ${params.microrna_seed_length} \\
        --weight_matrix ${params.weight_matrix}
    """
}