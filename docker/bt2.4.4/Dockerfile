FROM continuumio/miniconda3:4.10.3
MAINTAINER Minghao Chia <chiaminghao@gmail.com>

#update conda channel with bioconda and conda-forge. install bowtie2 and clean up to reduce size of docker image
RUN conda config --add channels conda-forge && conda config --add channels bioconda && \
conda install -c bioconda bowtie2=2.4.4 && \
conda install -c bioconda bbmap=38.93 && \
conda install -c bioconda samtools=1.13 && conda clean -afy
