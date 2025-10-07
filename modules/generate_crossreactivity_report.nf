process GENERATE_CROSSREACTIVITY_REPORT {
    tag "${params.run_id} - $gene_id - Generate Cross-Reactivity Report"
    // publishDir "${params.outdir}/${params.run_id}/${gene_id}", mode: 'copy'

    input:
    tuple val(gene_id), path(json_file)

    output:
    path "${gene_id}.crossreactivity.tsv", emit: crossreactivity_report

    script:
    def output_tsv = "${gene_id}.crossreactivity.tsv"
    """
    generate_crossreactivity_report.py \\
        --json ${json_file} \\
        --output ${output_tsv} \\
        --geneid_accession ${params.geneid_accession}
    """
}