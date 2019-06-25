# Base Image
FROM biocontainers/biocontainers:v1.0.0_cv4

# Metadata
LABEL base.image="biocontainers/biocontainers:v1.0.0_cv4"
LABEL version="1"
LABEL software="bcftools"
LABEL software.version="1.3.1"
LABEL description="Docker container for GWAS preprocessing"
LABEL website="https://samtools.github.io/bcftools/"
LABEL documentation="https://samtools.github.io/bcftools/"
LABEL license="https://samtools.github.io/bcftools/"
LABEL tags="Genomics"

USER root

RUN conda install bcftools=1.3.1 && \
    conda install htslib=1.9

RUN pip install pandas && \
    pip install numpy

ENV PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/opt/conda/bin:/home/biodocker/bin:/data/