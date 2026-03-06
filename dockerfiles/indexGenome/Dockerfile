# Use the continuumio/miniconda3 base image
FROM continuumio/miniconda3:23.3.1-0

RUN apt-get update && \
    apt-get install -y \
        bwa=0.7.17-6+b1 \
        samtools=1.11-1 && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda install -y \
        picard==2.18.7 && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* && \
    apt-get remove -y && \
    apt-get autoremove -y && \
    rm -rf /var/lib/apt/lists/* && \
    conda clean --all

# Ensure conda-installed binaries are on PATH for non-interactive shells
ENV PATH="/opt/conda/bin:$PATH"

# Metadata indicating the author or maintainer of the image
LABEL maintainer="Peter J Freeman <peter.j.freeman@manchester.ac.uk>"

# Metadata about the purpose or description of the image
LABEL description="Docker image containing bwa=0.7.17, bcftools=1.11, picard=2.18.7"

# Default command to run when the container starts
CMD ["/bin/bash"]
