process GENERATE_REPORT {
    tag "${params.run_id} - $gene_id - Generate Report"
    publishDir "${params.outdir}/${params.run_id}/${gene_id}", mode: 'copy'

    container 'python:3.10'

    input:
    tuple val(gene_id), path(json_file)

    output:
    path "${json_file.baseName}.tsv", emit: tsv

    script:
    def output_tsv = "${json_file.baseName}.tsv"
    """
    generate_report.py --json ${json_file} --output ${output_tsv}
    """
}