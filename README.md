# OLIGO-FINDER-NF: A Parallel RNAi Oligo Finder Pipeline

## Introduction

`OLIGO-FINDER-NF` is a scalable, parallel bioinformatics pipeline built with Nextflow. It is designed to identify potential RNAi oligo candidates from a set of target genes and screen them for off-target effects by aligning them against a comprehensive reference genome or transcriptome.

The pipeline takes a multi-FASTA file containing one or more target genes, generates candidate oligos based on user-defined criteria (length, GC content, etc.), aligns them using Bowtie, and produces a final, organized report for each gene, detailing the specificity of each oligo.

## Pipeline Workflow

The pipeline performs the following steps for each gene in the input file, executing them in parallel whenever possible:

1. **SPLIT_FASTA**: The input multi-FASTA file is split into individual FASTA files, one for each gene.

2. **GENERATE_OLIGO_CANDIDATE**: For each gene, a set of potential oligo sequences is generated based on specified length, GC content, sequence offsets, and forbidden motifs.

3. **BOWTIE_ALIGN**: The generated oligo candidates are aligned against a specified Bowtie index of a reference genome/transcriptome to find potential off-target matches.

4. **PARSE_SAM**: The raw SAM alignment output from Bowtie is parsed into a structured JSON format, aggregating matches by oligo and mismatch level.

5. **GENERATE_REPORT**: The final, human-readable TSV report is generated from the JSON file, summarizing the alignment results for each oligo candidate.

## Requirements

To run this pipeline, you will need:

1. **Nextflow**: [Installation instructions](https://www.nextflow.io/docs/latest/install.html#installation)

2. **A container engine**: 
    - **Docker**: The pipeline is configured to use Docker by default. [Installation instructions](https://docs.docker.com/get-started/get-docker/)
    - Singularity or Podman are also supported by Nextflow.

## Setup & Configuration

1. **Clone the repository**:

```bash
git clone <https://github.com/JulianH89/oligo-finder-nf.git>
cd OLIGO-FINDER-NF
```