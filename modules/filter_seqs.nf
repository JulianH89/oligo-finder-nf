process FILTER_SEQS {
    tag "${params.run_id} - $gene_id - Filter Sequences"
    publishDir "${params.outdir}/${params.run_id}/${gene_id}", mode: 'copy'

    input:
    tuple val(gene_id), path(seq)

    output:
    tuple val(gene_id), path("${gene_id}_filtered_seqs.tsv"), emit: filtered_seqs

    script:
    def output_seqs = "${gene_id}_filtered_seqs.tsv"

    """
    filter_sequences.py \
        --seq_file ${seq} \
        --min_gc ${params.min_gc} \
        --max_gc ${params.max_gc} \
        --microrna_hits_threshold ${params.microrna_hits_threshold} \
        --forbidden_motifs ${params.forbidden_motifs} \
        --output_file ${output_seqs}
    """
}