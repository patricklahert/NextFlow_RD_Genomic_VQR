# Use the continuumio/miniconda3 base image
FROM continuumio/miniconda3:23.3.1-0

RUN conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda install -y haplogrep==2.4.0 bcftools==1.21 plink2 && \
    conda clean --all -f -y

# Ensure conda-installed binaries are on PATH for non-interactive shells
ENV PATH="/opt/conda/bin:$PATH"

# Metadata indicating the author or maintainer of the image
LABEL maintainer="Peter J Freeman <peter.j.freeman@manchester.ac.uk>" \
      description="Docker image containing bcftools, haplogrep2, and plink2 for DNA identity analysis."

# Default command to run when the container starts
CMD ["/bin/bash", "-l"]
