FROM continuumio/miniconda3:4.10.3
MAINTAINER Minghao Chia <chiaminghao@gmail.com>

#update conda channel with bioconda and conda-forge. install software and clean up to reduce size of docker image
RUN conda config --add channels conda-forge && conda config --add channels bioconda && \
conda install -c bioconda eggnog-mapper=2.1.6 && \
conda install -c bioconda bbmap=38.93 && \
conda install -c bioconda bedtools=2.30.0 && conda clean -afy