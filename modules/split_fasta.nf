process SPLIT_FASTA {
    tag "${params.run_id} - Splitting Fasta File"

    container 'ubuntu:22.04'

    input:
    path target_gene

    output:
    path "*.fa", emit: fasta_files

    script:
    """
    awk -f ${baseDir}/bin/split_fasta.awk ${target_gene}
    """
}