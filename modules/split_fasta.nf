process SPLIT_FASTA {
    tag "${params.run_id} - Splitting Fasta File"
    publishDir "${params.outdir}/${params.run_id}", mode: 'copy'

    input:
    path split_script
    path target_gene

    output:
    path "split_fasta/*.fa", emit: fasta_files

    script:
    """
    mkdir -p split_fasta
    awk -v outdir="split_fasta" -f ${split_script} ${target_gene}
    """
}