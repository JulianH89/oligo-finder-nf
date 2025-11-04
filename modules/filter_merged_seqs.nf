process FILTER_MERGED_SEQS {
    tag "${params.run_id} - $gene_id - Filter Merged Sequences"
    container 'oligo-finder-env:latest'

    input:
    tuple val(gene_id), path(merged_seq)

    output:
    tuple val(gene_id), path("${gene_id}.filtered.tsv"), emit: filtered_seqs

    script:
    def output_seqs = "${gene_id}.filtered.tsv"

    """
    filter_sequences.py \
        --seq_file ${merged_seq} \
        --min_gc ${params.min_gc} \
        --max_gc ${params.max_gc} \
        --microrna_hits_threshold ${params.microrna_hits_threshold} \
        --forbidden_motifs ${params.forbidden_motifs} \
        --output_file ${output_seqs}
    """
}