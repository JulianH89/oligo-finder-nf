process MERGE_RESULTS {
    tag "${params.run_id} - $gene_id - Merge Results"

    input:
    tuple val(gene_id), path(filtered_metadata)
    path crossreactivity_report
    
    output:
    tuple val(gene_id), path("${gene_id}_merged.tsv"), emit: merged_result

    script:
    def output_tsv = "${gene_id}_merged.tsv"
    """
    merge_results.py \\
        --filtered_metadata ${filtered_metadata} \\
        --crossreactivity_report ${crossreactivity_report} \\
        --output ${output_tsv}
    """
}