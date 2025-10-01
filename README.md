# OLIGO-FINDER-NF: A Parallel RNAi Oligo Finder Pipeline

## Introduction

`OLIGO-FINDER-NF` is a scalable, parallel bioinformatics pipeline built with Nextflow. It is designed to identify potential RNAi oligo candidates from a set of target genes and screen them for off-target effects by aligning them against a comprehensive reference genome or transcriptome.

The pipeline takes a multi-FASTA file containing one or more target genes, generates candidate oligos based on user-defined criteria (length, GC content, etc.), aligns them using Bowtie, and produces a final, organized report for each gene, detailing the specificity of each oligo.

## Pipeline Workflow

The pipeline performs the following steps for each gene in the input file, executing them in parallel whenever possible:

1. **SPLIT_FASTA**: The input multi-FASTA file is split into individual FASTA files, one for each gene.

2. **GENERATE_METADATA**: For each gene, a set of sequences is generated, including surrounding sequence, oligos, refseq seeds, and reverse complement of oligos, etc.

3. **FILTER_METADATA**: Filter sequences based on GC content, microRNA hits and forbidden motifs.

4. **BOWTIE_ALIGN**: The generated oligo candidates are aligned against a specified Bowtie index of a reference genome/transcriptome to find potential off-target matches.

5. **PARSE_SAM**: The raw SAM alignment output from Bowtie is parsed into a structured JSON format, aggregating matches by oligo and mismatch level.

6. **GENERATE_CROSSREACTIVITY_REPORT**: The JOSN is used to generate the crossreacivity
 report of target genes.

7. **MERGE_RESULTS**: The filtered metadata file and crossreacivity report are merge to generate.

8. **CONVERT_TO_ORDER**: The merged reports are parsed to order format.

## Requirements

To run this pipeline, you will need:

1. **Nextflow**: [Installation instructions](https://www.nextflow.io/docs/latest/install.html#installation)

2. **A container engine**: 
    - **Docker**: The pipeline is configured to use Docker by default. [Installation instructions](https://docs.docker.com/get-started/get-docker/)

## Setup & Configuration

1. **Clone the repository**:

```bash
git clone <https://github.com/JulianH89/oligo-finder-nf.git>
cd OLIGO-FINDER-NF
```

2. **Prepare your data**:

    - **Target Genes**: FASTA file containing the genes of interest.

    - **Bowtie Index**: Pre-built Bowtie index files.

3. **Configure the pipeline**:

Open the `nextflow.config` file and edit the `params` block to match your file locations and desired settings.

### Core Parameters

| Parameter | Type | Default Value | Description |
|----------|----------|----------|----------|
| `run_id` | String |  | A unique name for the pipeline run. Used for organizing output. |
| `target_gene` | String(Path) |  | Path to the input multi-FASTA file containing target genes. |
| `outdir` | String(Path) | `$baseDir/results` | Path to the directory where all results and logs will be saved. |

### Reference Genome Parameters

| Parameter | Type | Default Value | Description |
|----------|----------|----------|----------|
| `bowtie_index_dir` | String(Path) |  | Path to the directory containing the Bowtie index files. |
| `bowtie_index_prefix` | String |  | The basename/prefix of the Bowtie index files (e.g., 'prefix' for `prefix.1.ebwt`). |

### Oligo Design Parameters

| Parameter | Type | Default Value | Description |
|----------|----------|----------|----------|
| `oligo_length` | Integer | `16` | The desired length of the oligo candidates to be generated. |
| `offset_5_prime` | Integer | `0` | Number of bases to trim from the 5' end of the target gene sequence before generating oligos. |
| `offset_3_prime` | Integer | `0` | Number of bases to trim from the 3' end of the target gene sequence. |
| `min_gc` | Float | `40.0` | The minimum allowed GC content percentage for an oligo candidate. |
| `max_gc` | Float | `60.0` | The maximum allowed GC content percentage for an oligo candidate. |
| `forbidden_motifs` | String | `GGG` | A comma-separated list of motifs that are not allowed in oligo candidates (e.g., `"GGG,AAAA"`). |

### Alignment Parameters

| Parameter | Type | Default Value | Description |
|----------|----------|----------|----------|
| `max_mismatch` | Integer | `3` | The maximum number of mismatches allowed during the Bowtie alignment (`-v` parameter). |

## Usage

To run the pipeline, execute the following command from the root directory of the project:

```bash
nextflow run main.nf
```

You can override any parameter from the command line using a double-dash prefix:

```bash
nextflow run main.nf --target_gene 'path/to/your/genes.fa' --run_id 'MyNewRun'
```

## Output

The pipeline will create an output directory specified by `params.outdir` (default is `results/`). The results are organized by run ID and then by gene ID.

### Directory Structure

```bash
results/
└── <run_id>/
    ├── <gene_A>/
    │   ├── gene_A.oligos.fa
    │   ├── gene_A.oligos.sam
    │   ├── gene_A.oligos.json
    │   └── gene_A.oligos.tsv  // Final Report for Gene A
    │
    ├── <gene_B>/
    │   ├── gene_B.oligos.fa
    │   ├── gene_B.oligos.sam
    │   └── ...
    │
    └── ...

```

### Final Report (`.tsv` file)

The final report is a tab-separated file with the following columns:

| Column | Description |
|----------|----------|
| Oligo_id | The unique identifier for the oligo candidate. |
| sequence | The DNA sequence of the oligo. |
| gc_content | The GC content percentage of the oligo. |
| mismatches | The number of mismatches (0, 1, 2, 3) for this alignment. The results will be displayed at field NM in SAM file. |
| matched_accession | A comma-separated list of reference sequences the oligo matched. |
| num_of_matched | The total count of reference sequences matched at this mismatch level. |

## Core Tools

This pipeline relies on the following core tools:

- **Nextflow**: Workflow management.

- **Bowtie**: Short read alignment.

- **Python**: Data processing scripts.

- **Docker**: Containerization.