# OLIGO-FINDER-NF: A Parallel RNAi Oligo Finder Pipeline

## Introduction

`OLIGO-FINDER-NF` is a scalable, parallel bioinformatics pipeline built with Nextflow. It is designed to identify potential RNAi oligo candidates from a set of target genes and screen them for off-target effects by aligning them against a comprehensive reference genome or transcriptome.

The pipeline takes a multi-FASTA file containing one or more target genes, generates candidate oligos based on user-defined criteria (length, GC content, etc.), aligns them using Bowtie, and produces a final, organized report for each gene, detailing the specificity of each oligo.

## Pipeline Workflow

The pipeline performs the following steps for each gene in the input file, executing them in parallel whenever possible:

1. **SPLIT_FASTA**: The input multi-FASTA file is split into individual FASTA files, one for each gene.

2. **GENERATE_SEQS**: For each gene, a set of sequences is generated, including surrounding sequence, oligos, refseq seeds, and reverse complement of oligos, etc.

3. **FILTER_SEQS**: Filter sequences based on GC content, microRNA hits and forbidden motifs.

4. **BOWTIE_ALIGN**: Align the generated refseq seeds against a reference genome/transcriptome to find off-target matches.

5. **PARSE_SAM**: Parse the SAM file for each gene into a structured JSON format.

6. **GENERATE_CROSSREACTIVITY_REPORT**: Generate the final TSV report for each gene from the JSON file.

7. **MERGE_RESULTS**: Merge the filtered sequences and cross-reactivity reports for each gene.

8. **CONVERT_TO_ORDER**: Convert the final report into a synthesis order format.

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

    - **Weight Matrix**: Weight Matrix.

    - **MicroRNA Seeds**: MicroRNA Seeds.

3. **Configure the pipeline**:

Open the `nextflow.config` file and edit the `params` block to match your file locations and desired settings.

### Parameters

#### Pipeline Parameters

| Parameter | Type | Default Value | Description |
|----------|----------|----------|----------|
| `run_id` | String |  | A unique name for the pipeline run. Used for organizing output. |
| `outdir` | String(Path) | `$baseDir/results` | Path to the directory where all results and logs will be saved. |

#### Target Gene Parameters

| Parameter | Type | Default Value | Description |
|----------|----------|----------|----------|
| `target_gene` | String(Path) |  | Path to the input multi-FASTA file containing target genes. |

#### Reference Genome Parameters

| Parameter | Type | Default Value | Description |
|----------|----------|----------|----------|
| `bowtie_index_dir` | String(Path) |  | Path to the directory containing the Bowtie index files. |
| `bowtie_index_prefix` | String |  | The basename/prefix of the Bowtie index files (e.g., 'prefix' for `prefix.1.ebwt`). |

#### Files Parameters

| Parameter | Type | Default Value | Description |
|----------|----------|----------|----------|
| `weight_matrix` | String(Path) |  | Path to the Weight Matrix files. |
| `microrna_seeds` | String(Path) |  | Path to the microRNA Seeds files. |

#### sdRNAi Design Parameters

| Parameter | Type | Default Value | Description |
|----------|----------|----------|----------|
| `surrounding_region_length` | Integer | `45` | The desired length of the surrounding regions |
| `oligo_length` | Integer | `20` | The desired length of the oligo candidates to be generated. |
| `offset_5_prime` | Integer | `16` | Number of bases to trim from the 5' start of the surrounding regions. |
| `offset_refseq_seed` | Integer | `3` | Number of bases to trim from the 5' start of the oligos. |
| `refseq_seed_length` | Integer | `16` | The desired length of the refseq seeds to be generated. |
| `offset_microrna` | Integer | `3` | Number of bases to trim from the 5' start of the microRNA seeds. |
| `microrna_seed_length` | Integer | `7` | The desired length of the microRNA seeds to be generated. |

#### Filtering parameters

| Parameter | Type | Default Value | Description |
|----------|----------|----------|----------|
| `min_gc` | Float | `40.0` | The minimum allowed GC content percentage for an oligo candidate. |
| `max_gc` | Float | `60.0` | The maximum allowed GC content percentage for an oligo candidate. |
| `microrna_hits_threshold` | String | `1` | The maximum allowed microRNA hits for am oligo candidate. |
| `forbidden_motifs` | String | `GGG` | A comma-separated list of motifs that are not allowed in oligo candidates (e.g., `"GGG,AAAA"`). |

#### Alignment Parameters

| Parameter | Type | Default Value | Description |
|----------|----------|----------|----------|
| `max_mismatch` | Integer | `3` | The maximum number of mismatches allowed during the Bowtie alignment (`-v` parameter). |

#### Synthesis Order Parameters

| Parameter | Type | Default Value | Description |
|----------|----------|----------|----------|
| `sense_length` | Integer | `14` | The desired length of the sense length. |
| `antisense_length` | Integer | `19` | The desired length of the antisense length. |

## Usage

To run the pipeline, execute the following command from the root directory of the project:

```bash
nextflow run main.nf -profile docker
```

You can override any parameter from the command line using a double-dash prefix:

```bash
nextflow run main.nf -profile docker --run_id 'My_New_Run'
```

## Output

The pipeline will create an output directory specified by `params.outdir` (default is `results/`). The results are organized by run ID and then by gene ID.

### Directory Structure

```bash
results/
└── <run_id>/
    ├── <gene_A>/
    │   ├── gene_A_filtered_seqs.tsv
    │   ├── gene_A_report.tsv
    │   ├── gene_A.crossreactivity.tsv
    │   ├── gene_A.json
    │   ├── gene_A.order.tsv
    │   └── gene.seqs.tsv
    │
    ├── <gene_B>/
    │   ├── gene_B_filtered_seqs.tsv
    │   ├── gene_B_report.tsv
    │   └── ...
    │
    └── ...

```

#### Explanation of output files

| File name | Description |
|----------|----------|
| `*_filtered_seqs` |  |
| `*_report.tsv` |  |
| `*.crossreactivity.tsv` |  |
| `*.json` |  |
| `*.order.tsv` |  |
| `*.seqs.tsv` |  |

### Final Report (`.order.tsv` file)

The final report is a tab-separated file with the following columns:

| Column | Description |
|----------|----------|
| ID | The unique identifier for the oligo candidate. |
| Surrounding_Region | The DNA sequence of he surrounding region. |
| Oligo | The DNA sequence of the oligo. |
| Oligo_RC | The DNA sequence of the reverse complement of the oligo. |
| GC_Content | The GC content percentage of the oligo. |
| Sense_Tripurine |  |
| Antisense_Tripurine |  |
| Sense_FM |  |
| Antisense_FM |  |
| Refseq_Seed | The DNA sequence of the refseq seed(the sequence actually map to the reference). |
| Score | The score is calculated based on surrounding region and weight matrix. |
| mismatch_level | The number of mismatches (0, 1, 2, 3) for this alignment. The results will be displayed at field NM in SAM file. |
| num_of_matched_accessions | The total count of reference sequences matched at this mismatch level. |
| MicroRNA_Seed | The DNA sequence of microRNA seed of the oligo. |
| MicroRNA_Hits | The hits of the microRNA seed againt database. |
| matched_accession | A comma-separated list of reference sequences the oligo matched. |

## Core Tools

This pipeline relies on the following core tools:

- **Nextflow**: Workflow management.

- **Bowtie**: Short read alignment.

- **Python**: Data processing scripts.

- **Docker**: Containerization.