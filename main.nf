#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
========================================================================================
    OLIGO-FINDER-NF PIPELINE
========================================================================================
    Pipeline for aligning potential RNAi oligos against a reference gene set using Bowtie.
----------------------------------------------------------------------------------------
*/

// --- HELPER FUNCTION ---
// Function to validate required parameters.
def validate_params() {
    if (!params.run_id) {
        error "ERROR: A run ID must be provided using --run_id <ID>"
    }
    if (!params.oligos) {
        error "ERROR: A query file must be provided using --oligos <path/to/file.fa>"
    }
    if (!params.bowtie_index_dir) {
        error "ERROR: A Bowtie index directory must be provided using --bowtie_index_dir <path/to/dir>"
    }
    if (!params.bowtie_index_prefix) {
        error "ERROR: A Bowtie index prefix must be provided using --bowtie_index_prefix <prefix>"
    }
    if (!params.max_mismatch) {
        error "ERROR: A maximum number of mismatches must be provided using --max_mismatch <int>"
    }
}

// --- LOG ---
log.info """
         O L I G O T I E - RNAi Oligo Finder
         ===================================
         Run ID         : ${params.run_id}
         Reference      : ${params.bowtie_index_dir}/${params.bowtie_index_prefix}
         Oligo File     : ${params.oligos}
         Oligo Length   : ${params.'oligo_length'}
         Max Mismatches : ${params.max_mismatch}
         Output Dir     : ${params.outdir}
         """
         .stripIndent()


// --- MODULES ---
include { BOWTIE_ALIGN } from './modules/bowtie_align'
include { PARSE_SAM }    from './modules/parse_sam'
include { GENERATE_REPORT }   from './modules/generate_report'


// --- WORKFLOW ---
workflow {

    // Validate parameters at the start of the workflow.
    validate_params()

    // Create a channel that emits the path to the Bowtie index.
    Channel
        .of([ file(params.bowtie_index_dir, checkIfExists: true), params.bowtie_index_prefix ])
        .set { ch_bowtie_index }

    // Create a channel for the oligo sequences.
    Channel
        .fromPath(params.oligos, checkIfExists: true)
        .set { ch_oligos }


    // 1. Align the oligo sequences.
    BOWTIE_ALIGN (
        params.run_id,
        ch_oligos,
        ch_bowtie_index,
        params.max_mismatch
    )

    // 2. Parse the SAM file into a structured JSON format
    PARSE_SAM (
        params.run_id,
        BOWTIE_ALIGN.out.sam
    )

    // 3. Generate the final TSV report from the JSON file
    GENERATE_REPORT (
        params.run_id,
        PARSE_SAM.out.json
    )

    // Log the final report path for easy viewing.
    GENERATE_REPORT.out.tsv.view()
}