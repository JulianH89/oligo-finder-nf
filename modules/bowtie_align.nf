process BOWTIE_ALIGN {
    tag "${params.run_id} - $gene_id - Bowtie Alignment"
    publishDir "${params.outdir}/${params.run_id}/${gene_id}", mode: 'copy'

    input:
    tuple val(gene_id), path(metadata_seq)

    output:
    tuple val(gene_id), path("${gene_id}.sam"), emit: sam

    script:
    def threads = task.cpus
    def refseq_seed_fasta = "${gene_id}_refseq_seed.fasta"
    def output_sam = "${gene_id}.sam"
    def bowtie_index_path = "${params.bowtie_index_dir}/${params.bowtie_index_prefix}"

    // Bowtie command to perform the alignment.
    // Allowing up to 'max_mismatch' mismatches.
    // The --norc option is used to prevent alignment to the reverse complement strand.
    """
    # Extract oligo sequences from metadata file
    awk 'NR>1 {print ">"\$1"\\n"\$5}' ${metadata_seq} > ${refseq_seed_fasta}
    
    bowtie --threads ${threads} --quiet -a --norc \\
        ${bowtie_index_path} \\
        -f ${refseq_seed_fasta} \\
        -S ${output_sam} \\
        -v ${params.max_mismatch}
    """
}