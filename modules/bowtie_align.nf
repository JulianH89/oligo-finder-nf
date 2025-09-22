process BOWTIE_ALIGN {
    tag "${params.run_id} - $gene_id - Bowtie Alignment"
    publishDir "${params.outdir}/${params.run_id}/${gene_id}", mode: 'copy'

    container 'biocontainers/bowtie:v1.1.2_cv4'

    input:
    tuple val(gene_id), path(oligos_fasta)

    output:
    tuple val(gene_id), path("${oligos_fasta.baseName}.sam"), emit: sam

    script:
    def threads = task.cpus
    def output_sam = "${oligos_fasta.baseName}.sam"
    def bowtie_index_path = "${params.bowtie_index_dir}/${params.bowtie_index_prefix}"

    // Bowtie command to perform the alignment.
    // Allowing up to 'max_mismatch' mismatches.
    // The --norc option is used to prevent alignment to the reverse complement strand.
    """
    bowtie --threads ${threads} --quiet -a --norc \\
        ${bowtie_index_path} \\
        -f ${oligos_fasta} \\
        -S ${output_sam} \\
        -v ${params.max_mismatch}
    """
}