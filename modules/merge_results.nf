process MERGE_RESULTS {
    tag "${params.run_id} - $gene_id - Merge Results"
    container 'oligo-finder-env:latest'

    input:
    tuple val(gene_id), path(metadata)
    path crossreactivity_report
    
    output:
    tuple val(gene_id), path("${gene_id}.compete.tsv"), emit: merged_result

    script:
    def output_tsv = "${gene_id}.compete.tsv"
    """
    merge_results.py \\
        --filtered_metadata ${metadata} \\
        --crossreactivity_report ${crossreactivity_report} \\
        --output ${output_tsv}
    """
}