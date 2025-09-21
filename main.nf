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
    if (!params.target_gene) {
        error "ERROR: A target gene FASTA must be provided using --target_gene <path/to/file.fa>"
    }
    if (!params.bowtie_index_dir) {
        error "ERROR: A Bowtie index directory must be provided using --bowtie_index_dir <path/to/dir>"
    }
    if (!params.bowtie_index_prefix) {
        error "ERROR: A Bowtie index prefix must be provided using --bowtie_index_prefix <prefix>"
    }
    if (!params.oligo_length) {
        error "ERROR: An oligo length must be provided using --oligo_length <int>"
    }
}

// --- LOG ---
log.info """
         O L I G O T I E - RNAi Oligo Finder
         ===================================
         Run ID             : ${params.run_id}
         Reference          : ${params.bowtie_index_dir}/${params.bowtie_index_prefix}
         Target Gene        : ${params.target_gene_fasta}
         Oligo Length       : ${params.oligo_length}
         5' Offset          : ${params.offset_5_prime}
         3' Offset          : ${params.offset_3_prime}
         GC Range           : ${params.min_gc}% - ${params.max_gc}%
         Forbidden Motifs   : ${params.forbidden_motifs}
         Max Mismatches     : ${params.max_mismatch}
         Output Dir         : ${params.outdir}
         """
         .stripIndent()


// --- MODULES ---
include { GENERATE_OLIGO_CANDIDATE } from './modules/generate_oligo_candidate'
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

    // 0. Generate oligo candidates from the target gene
    GENERATE_OLIGO_CANDIDATE (
        params.run_id,
        file(params.target_gene, checkIfExists: true),
        params.oligo_length,
        params.offset_5_prime,
        params.offset_3_prime,
        params.min_gc,
        params.max_gc,
        params.forbidden_motifs
    )

    // 1. Align the oligo sequences.
    BOWTIE_ALIGN (
        params.run_id,
        ch_bowtie_index,
        GENERATE_OLIGO_CANDIDATE.out.oligos_fasta,
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

}

