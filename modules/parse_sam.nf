process PARSE_SAM {
    tag "$run_id"
    publishDir "${params.outdir}/parsed_result/${run_id}", mode: 'copy'

    container 'python:3.10-slim'

    input:
    val run_id
    path sam_file

    output:
    path "${sam_file.baseName}.json", emit: json

    script:
    def output_json = "${sam_file.baseName}.json"
    """
    parse_sam.py --sam ${sam_file} --output ${output_json}
    """
}