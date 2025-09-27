process FILTER_METADATA {
    tag "${params.run_id} - $gene_id - Filter Metadata"
    publishDir "${params.outdir}/${params.run_id}/${gene_id}", mode: 'copy'

    input:
    tuple val(gene_id), path(seq_metadata)

    output:
    tuple val(gene_id), path("${gene_id}_filtered_metadata.tsv"), emit: filtered_metadata

    script:
    def output_metadata = "${gene_id}_filtered_metadata.tsv"

    """
    filter_metadata.py \
        --metadata_file ${seq_metadata} \
        --min_gc ${params.min_gc} \
        --max_gc ${params.max_gc} \
        --microrna_hits_threshold ${params.microrna_hits_threshold} \
        --forbidden_motifs ${params.forbidden_motifs} \
        --output_file ${output_metadata}
    """
}