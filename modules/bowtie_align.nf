process BOWTIE_ALIGN {
    tag "$run_id"
    publishDir "results/alignment", mode: 'copy'

    container 'biocontainers/bowtie:v1.2.2dfsg-4-deb_cv1'

    input:
    tuple val(run_id), path(query_seq)
    path index_dir
    val index_prefix
    val max_mismatch                    // The value for the -v parameter

    output:
    path "${run_id}.sam", emit: sam  // The output SAM file, emitted to a channel named 'sam'

    script:
    def threads = task.cpus
    def output_sam = "${run_id}.sam"
    def bowtie_index_path = "${index_dir}/${index_prefix}"

    """
    bowtie --threads $threads --quiet -a \
        ${bowtie_index_path} \
        -f ${query_seq} \
        -S ${output_sam} \
        -v ${max_mismatch}
    """
}