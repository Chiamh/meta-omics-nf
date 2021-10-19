FROM continuumio/miniconda3:4.10.3
MAINTAINER Minghao Chia <chiaminghao@gmail.com>

#update conda channel with bioconda and conda-forge. install software and clean up to reduce size of docker image
RUN conda config --add channels conda-forge && conda config --add channels bioconda && \
conda install -c bioconda kraken2=2.1.2 && \
conda install -c bioconda bracken=2.6.1 && \ 
conda install -c bioconda diamond=2.0.12 && \
conda install -c bioconda star=2.7.9a && \
conda install -c bioconda fastp=0.22.0 && \
conda install -c bioconda bbmap=38.93 && \
conda install -c bioconda bwa=0.7.17 && conda clean -afy 
