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
        Refseq Seed Offset         : ${params.offset_refseq_seed}
        Refseq Seed Length         : ${params.refseq_seed_length}
        MicroRNA Seeds             : ${params.microrna_seeds}
        MicroRNA Offset            : ${params.offset_microrna}
        MicroRNA Length            : ${params.microrna_seed_length}
        GC Range                   : ${params.min_gc}% - ${params.max_gc}%
        Forbidden Motifs           : ${params.forbidden_motifs}
        Max Mismatches             : ${params.max_mismatch}
        Sense Length               : ${params.sense_length}
        Antisense Length           : ${params.antisense_length}
        Output Dir                 : ${params.outdir}/${params.run_id}/
        ===================================
        """
        .stripIndent()


// --- MODULES ---
include { SPLIT_FASTA } from './modules/split_fasta'
include { GENERATE_SEQS } from './modules/generate_seqs'
include { FILTER_SEQS } from './modules/filter_seqs'
include { BOWTIE_ALIGN } from './modules/bowtie_align'
include { PARSE_SAM }    from './modules/parse_sam'
include { GENERATE_CROSSREACTIVITY_REPORT }   from './modules/generate_crossreactivity_report'
include { MERGE_RESULTS } from './modules/merge_results'
include { CONVERT_TO_ORDER } from './modules/convert_to_order'



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
    GENERATE_SEQS (
        ch_genes
    )

    // 2. Filter the metadata for each gene in parallel
    FILTER_SEQS (
        GENERATE_SEQS.out.seq
    )

    // 3. Align the oligo sequences for each gene.
    BOWTIE_ALIGN (
        FILTER_SEQS.out.filtered_seqs
    )

    // 4. Parse the SAM file for each gene into a structured JSON format
    PARSE_SAM (
        BOWTIE_ALIGN.out.sam
    )

    // 5. Generate the final TSV report for each gene from the JSON file
    GENERATE_CROSSREACTIVITY_REPORT (
        PARSE_SAM.out.json
    )

    // 6. Merge the filtered sequences and cross-reactivity reports for each gene
    MERGE_RESULTS (
        FILTER_SEQS.out.filtered_seqs,
        GENERATE_CROSSREACTIVITY_REPORT.out.crossreactivity_report
    )

    // 7. Convert the final report into a synthesis order format
    CONVERT_TO_ORDER (
        MERGE_RESULTS.out.merged_result
    )

}

