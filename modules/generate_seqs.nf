process GENERATE_SEQS {
    tag "${params.run_id} - $gene_id - Generate Sequences"
    publishDir "${params.outdir}/${params.run_id}/${gene_id}", mode: 'copy'

    input:
    tuple val(gene_id), path(target_gene)

    output:
    tuple val(gene_id), path("${gene_id}.seqs.tsv"), emit: seq

    script:
    def seq = "${gene_id}.seqs.tsv"
    """
    generate_sequences.py \\
        --input_fasta ${target_gene} \\
        --output ${seq} \\
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