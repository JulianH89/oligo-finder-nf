process GENERATE_REPORT {
    tag "$run_id - Generate Report"
    publishDir "${params.outdir}/${run_id}", mode: 'copy'

    container 'python:3.10-slim'

    input:
    val run_id
    path json_file

    output:
    path "${json_file.baseName}.tsv", emit: tsv

    script:
    def output_tsv = "${json_file.baseName}.tsv"
    """
    generate_report.py --json ${json_file} --output ${output_tsv}
    """
}