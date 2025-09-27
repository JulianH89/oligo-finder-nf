process MERGE_RESULTS {
    tag "${params.run_id} - $gene_id - Merge Results"
    publishDir "${params.outdir}/${params.run_id}/${gene_id}", mode: 'copy'

    input:
    tuple val(gene_id), path(filtered_metadata)
    path crossreactivity_report
    
    output:
    path "${gene_id}_report.tsv", emit: report

    script:
    def output_tsv = "${gene_id}_report.tsv"
    """
    merge_results.py \\
        --filtered_metadata ${filtered_metadata} \\
        --crossreactivity_report ${crossreactivity_report} \\
        --output ${output_tsv}
    """
}