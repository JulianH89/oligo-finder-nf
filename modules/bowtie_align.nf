process BOWTIE_ALIGN {
    tag "$run_id - Bowtie Alignment"
    publishDir "${params.outdir}/${run_id}", mode: 'copy'

    container 'biocontainers/bowtie:v1.1.2_cv4'

    input:
    val run_id
    tuple path(index_dir), val(index_prefix)
    path(oligos_fasta)
    val max_mismatch                    // The value for the -v parameter

    output:
    path "${oligos_fasta.baseName}.sam", emit: sam  // The output SAM file, emitted to a channel named 'sam'

    script:
    def threads = task.cpus
    def output_sam = "${oligos_fasta.baseName}.sam"
    def bowtie_index_path = "${index_dir}/${index_prefix}"

    // Bowtie command to perform the alignment.
    // Allowing up to 'max_mismatch' mismatches.
    // The --norc option is used to prevent alignment to the reverse complement strand.
    """
    bowtie --threads ${threads} --quiet -a --norc \
        ${bowtie_index_path} \
        -f ${oligos_fasta} \
        -S ${output_sam} \
        -v ${max_mismatch}
    """
}