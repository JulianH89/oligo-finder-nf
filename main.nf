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

// --- MODULES ---
include { SPLIT_FASTA } from './modules/split_fasta'
include { GENERATE_SEQS } from './modules/generate_seqs'
include { FILTER_SEQS } from './modules/filter_seqs'
include { BOWTIE_ALIGN } from './modules/bowtie_align'
include { PARSE_SAM }    from './modules/parse_sam'
include { GENERATE_CROSSREACTIVITY_REPORT }   from './modules/generate_crossreactivity_report'
include { MERGE_RESULTS } from './modules/merge_results'
include { GENERATE_FINAL_REPORT } from './modules/generate_final_report'



// --- WORKFLOW ---
workflow {

    // Validate parameters at the start of the workflow.
    validate_params()

    // ===== RECORD METADATA =====
    def metadata = [
        run_id: params.run_id,
        timestamp: new Date().format('yyyy-MM-dd HH:mm:ss'),
        nextflow_version: nextflow.version.toString(),
        workflow_session: workflow.sessionId,
        profile: workflow.profile,
    
        // Reference genome parameters
        bowtie_index_dir: params.bowtie_index_dir,
        bowtie_index_prefix: params.bowtie_index_prefix,
    
        // Input files
        target_gene: params.target_gene,
        weight_matrix: params.weight_matrix,
        microrna_seeds: params.microrna_seeds,
        geneid_accession: params.geneid_accession,
        cds_region: params.cds_region,
    
        // Design parameters
        surrounding_region_length: params.surrounding_region_length,
        oligo_length: params.oligo_length,
        offset_5_prime: params.offset_5_prime,
        offset_refseq_seed: params.offset_refseq_seed,
        refseq_seed_length: params.refseq_seed_length,
        offset_microrna: params.offset_microrna,
        microrna_seed_length: params.microrna_seed_length,
    
        // Filtering parameters
        min_gc: params.min_gc,
        max_gc: params.max_gc,
        microrna_hits_threshold: params.microrna_hits_threshold,
        forbidden_motifs: params.forbidden_motifs,
    
        // Alignment parameters
        max_mismatch: params.max_mismatch,
    
        // Synthesis parameters
        sense_length: params.sense_length,
        antisense_length: params.antisense_length,
    
        // Output directory
        outdir: params.outdir
    ]

    // Create output directory if it doesn't exist
    new File("${params.outdir}/${params.run_id}").mkdirs()

    // Write metadata to JSON file
    def metadata_json = file("${params.outdir}/${params.run_id}/run_metadata.json")
    metadata_json.text = groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(metadata))

    // ===== END: RECORD METADATA =====

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

    // 7. Generate the final report in a chemically-modified format
    GENERATE_FINAL_REPORT (
        MERGE_RESULTS.out.merged_result
    )

}

