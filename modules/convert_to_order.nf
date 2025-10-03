process CONVERT_TO_ORDER {
    tag "${params.run_id} - $gene_id - Generate chemically-modified format"
    publishDir "${params.outdir}/${params.run_id}/${gene_id}", mode: 'copy'

    input:
    tuple val(gene_id), path(report)

    output:
    path "${gene_id}.final.tsv", emit: final_report

    script:
    def output_order = "${gene_id}.final.tsv"
    """
    convert_to_order.py \\
        --report_tsv ${report} \\
        --output_tsv ${output_order} \\
        --sense_length ${params.sense_length} \\
        --antisense_length ${params.antisense_length}
    """
}