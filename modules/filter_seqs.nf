process FILTER_SEQS {
    tag "${params.run_id} - $gene_id - Filter Sequences"
    container 'oligo-finder-env:latest'

    input:
    tuple val(gene_id), path(seq)

    output:
    tuple val(gene_id), path("${gene_id}.filtered_seqs.tsv"), emit: filtered_seqs

    script:
    def output_seqs = "${gene_id}.filtered_seqs.tsv"

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