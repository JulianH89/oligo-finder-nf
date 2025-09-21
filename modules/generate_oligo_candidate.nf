process GENERATE_OLIGO_CANDIDATE {
    tag "$run_id - Generate Oligo Candidates"
    publishDir "${params.outdir}/${run_id}", mode: 'copy'

    container 'python:3.10-slim'

    input:
    val run_id
    path target_gene
    val oligo_length
    val offset_5
    val offset_3

    output:
    path "${target_gene.baseName}.oligos.fa", emit: oligos_fasta

    script:
    def output_fasta = "${target_gene.baseName}.oligos.fa"
    """
    generate_oligo_candidate.py \\
        --input-fasta ${target_gene} \\
        --output-fasta ${output_fasta} \\
        --oligo-length ${oligo_length} \\
        --offset-5 ${offset_5} \\
        --offset-3 ${offset_3}
    """
}