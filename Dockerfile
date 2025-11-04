FROM condaforge/mambaforge:latest

RUN mamba install --yes -c bioconda -c conda-forge \
    python=3.10 \
    pandas=2.1.0 \
    bowtie=1.3.1 \
    openpyxl