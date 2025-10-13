process SPLIT_FASTA {
    tag "${params.run_id} - Splitting Fasta File"

    container 'ubuntu:22.04'
    containerOptions '-u $(id -u):$(id -g)'

    input:
    path target_gene

    output:
    path "split_fasta/*.fa", emit: fasta_files

    script:
    """
    mkdir -p split_fasta
    awk -v outdir="split_fasta" -f ${baseDir}/bin/split_fasta.awk ${target_gene}
    """
}