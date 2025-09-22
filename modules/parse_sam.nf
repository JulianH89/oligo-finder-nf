process PARSE_SAM {
    tag "${params.run_id} - $gene_id - Parse SAM File"
    publishDir "${params.outdir}/${params.run_id}/${gene_id}", mode: 'copy'

    container 'python:3.10'

    input:
    tuple val(gene_id), path(sam_file)

    output:
    tuple val(gene_id), path("${sam_file.baseName}.json"), emit: json

    script:
    def output_json = "${sam_file.baseName}.json"
    """
    parse_sam.py --sam ${sam_file} --output ${output_json}
    """
}