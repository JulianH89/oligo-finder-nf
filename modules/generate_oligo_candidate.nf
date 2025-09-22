process GENERATE_OLIGO_CANDIDATE {
    tag "${params.run_id} - $gene_id - Generate Oligo Candidates"
    publishDir "${params.outdir}/${params.run_id}/${gene_id}", mode: 'copy'

    container 'python:3.10'

    input:
    tuple val(gene_id), path(target_gene)

    output:
    tuple val(gene_id), path("${target_gene.baseName}.oligos.fa"), emit: oligos_fasta

    script:
    def output_fasta = "${target_gene.baseName}.oligos.fa"
    """
    generate_oligo_candidate.py \\
        --input-fasta ${target_gene} \\
        --output-fasta ${output_fasta} \\
        --oligo-length ${params.oligo_length} \\
        --offset-5 ${params.offset_5_prime} \\
        --offset-3 ${params.offset_3_prime} \\
        --min-gc ${params.min_gc} \\
        --max-gc ${params.max_gc} \\
        --forbidden-motifs "${params.forbidden_motifs}"
    """
}