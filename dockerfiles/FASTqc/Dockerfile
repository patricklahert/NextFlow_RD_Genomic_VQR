# Use the continuumio/miniconda3 base image
FROM continuumio/miniconda3:23.3.1-0

RUN conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda install -y \
        fastqc==0.12.1 && \
    conda clean --all -f -y

# Ensure conda-installed binaries are on PATH for non-interactive shells
ENV PATH="/opt/conda/bin:$PATH"

# Metadata indicating the author or maintainer of the image
LABEL maintainer="Peter J Freeman <peter.j.freeman@manchester.ac.uk>"

# Metadata about the purpose or description of the image
LABEL description="Docker image containing FASTqc v0.12.1"

# Default command to run when the container starts
CMD ["/bin/bash", "-l"]