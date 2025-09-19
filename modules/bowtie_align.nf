process BOWTIE_ALIGN {
    tag "$run_id"
    publishDir "${params.outdir}/alignment/${run_id}", mode: 'copy'

    container 'biocontainers/bowtie:v1.2.2dfsg-4-deb_cv1'

    input:
    val run_id                          // The run ID for tagging and output directory
    path oligo
    tuple path(index_dir), val(index_prefix)
    val max_mismatch                    // The value for the -v parameter

    output:
    path "${oligo.baseName}.sam", emit: sam  // The output SAM file, emitted to a channel named 'sam'

    script:
    def threads = task.cpus
    def output_sam = "${oligo.baseName}.sam"
    def bowtie_index_path = "${index_dir}/${index_prefix}"

    """
    bowtie --threads $threads --quiet -a \
        ${bowtie_index_path} \
        -f ${oligo} \
        -S ${output_sam} \
        -v ${max_mismatch}
    """
}