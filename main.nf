#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
========================================================================================
    OLIGO-FINDER-NF PIPELINE
========================================================================================
    Pipeline for aligning potential RNAi oligos against a reference gene set using Bowtie.
----------------------------------------------------------------------------------------
*/

// --- LOG ---
log.info """
         O L I G O T I E - RNAi Oligo Finder
         ===================================
         Reference Dir  : ${params.bowtie_index_dir}
         Index Prefix   : ${params.bowtie_index_prefix}
         Oligo Files    : ${params.oligos}
         Mismatches     : ${params.max_mismatch}
         Output Dir     : ${params.outdir}
         """
         .stripIndent()


// --- MODULES ---
include { BOWTIE_ALIGN } from './modules/bowtie_align'


// --- WORKFLOW ---
workflow {

    // Create a channel that emits the path to the Bowtie index.
    Channel
        .fromPath(params.bowtie_index_dir, checkIfExists: true)
        .set { ch_bowtie_index_dir }

    // Create a channel for the oligo sequences.
    Channel
        .fromFilePairs(params.oligos)
        .ifEmpty { error "Cannot find any fasta files matching: ${params.oligos}" }
        .set { ch_oligos }


    // Call the BOWTIE_ALIGN module.
    // The inputs are provided in the order they are defined in the module's 'input' block.
    BOWTIE_ALIGN (
        ch_oligos,
        ch_bowtie_index_dir,
        params.bowtie_index_prefix,
        params.max_mismatch
    )

    // Log the output of the alignment process for easy viewing.
    BOWTIE_ALIGN.out.sam.view()
}