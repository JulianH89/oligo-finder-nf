#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
==========================================================================================
    OLIGO-FINDER-NF PIPELINE
==========================================================================================
    Pipeline for aligning potential RNAi oligos against a reference gene set using Bowtie.
------------------------------------------------------------------------------------------
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
        Run ID                     : ${params.run_id}
        Reference                  : ${params.bowtie_index_dir}/${params.bowtie_index_prefix}
        Target Gene                : ${params.target_gene}
        Surrounding Region Length  : ${params.surrounding_region_length}
        Oligo Length               : ${params.oligo_length}
        5' Offset                  : ${params.offset_5_prime}
        Gene Region Offset         : ${params.offset_gene_region}
        Gene Region Length         : ${params.gene_region_length}
        MicroRNA Offset            : ${params.offset_microrna}
        MicroRNA Length            : ${params.microrna_seed_length}
        GC Range                   : ${params.min_gc}% - ${params.max_gc}%
        Forbidden Motifs           : ${params.forbidden_motifs}
        Max Mismatches             : ${params.max_mismatch}
        Output Dir                 : ${params.outdir}/${params.run_id}/
        ===================================
        """
        .stripIndent()


// --- MODULES ---
include { SPLIT_FASTA } from './modules/split_fasta'
include { GENERATE_METADATA } from './modules/generate_metadata'
include { BOWTIE_ALIGN } from './modules/bowtie_align'
include { PARSE_SAM }    from './modules/parse_sam'
include { GENERATE_REPORT }   from './modules/generate_report'



// --- WORKFLOW ---
workflow {

    // Validate parameters at the start of the workflow.
    validate_params()

    // Create a value tuple for the Bowtie index
    def bowtie_index_tuple = tuple (
        file(params.bowtie_index_dir, checkIfExists: true), 
        params.bowtie_index_prefix
    )

    // 0. Split the multi-fasta file into a channel of single-gene fasta files
    SPLIT_FASTA (
        file(params.target_gene, checkIfExists: true)
    )

    // Take the channel of files from SPLIT_FASTA and transform it into the tuple format.
    SPLIT_FASTA.out.fasta_files
        .flatten()
        .map { file -> tuple(file.baseName, file) }
        .set { ch_genes }

    // 1. Generate metadata and oligo candidates from each target gene in parallel
    GENERATE_METADATA (
        ch_genes
    )

    // 2. Align the oligo sequences for each gene.
    BOWTIE_ALIGN (
        GENERATE_METADATA.out.seq_metadata
    )

    // 3. Parse the SAM file for each gene into a structured JSON format
    PARSE_SAM (
        BOWTIE_ALIGN.out.sam
    )

    // 4. Generate the final TSV report for each gene from the JSON file
    GENERATE_REPORT (
        PARSE_SAM.out.json
    )

}

