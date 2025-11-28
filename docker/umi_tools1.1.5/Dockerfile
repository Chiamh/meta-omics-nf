FROM continuumio/miniconda3:23.10.0-1

RUN conda config --add channels conda-forge && conda config --add channels bioconda && \
conda install -y -c bioconda umi_tools=1.1.5 && conda clean -afy
