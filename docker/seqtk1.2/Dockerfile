FROM continuumio/miniconda3:4.10.3
MAINTAINER Minghao Chia <chiaminghao@gmail.com>

#update conda channel with bioconda and conda-forge. install bowtie2 and clean up to reduce size of docker image
RUN conda config --add channels conda-forge && conda config --add channels bioconda && \
conda install -y -c bioconda samtools=1.13 && \
conda install -c bioconda seqtk=1.2 && conda clean -afy