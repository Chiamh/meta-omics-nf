FROM continuumio/miniconda3:4.10.3
MAINTAINER Minghao Chia <chiaminghao@gmail.com>

#update conda channel with bioconda and conda-forge. install software and clean up to reduce size of docker image
RUN conda config --add channels conda-forge && conda config --add channels bioconda && \
conda install -c bioconda vsearch=2.30.1 && conda clean -afy
