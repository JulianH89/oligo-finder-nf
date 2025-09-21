process PARSE_SAM {
    tag "$run_id - Parse SAM File"
    publishDir "${params.outdir}/${run_id}", mode: 'copy'

    container 'python:3.10'

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