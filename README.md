# Shotgun Metagenomics & Metatranscriptomics Nextflow Pipeline

## Introduction

Chiamh/meta-omics-nf is a bioinformatics pipeline that takes raw metagenomic and/or metatranscriptomic reads and annotates them at the taxonomic and functional levels.

The pipeline is built using [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) , a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. 
It uses Docker containers (also compatible with Singularity) for ease of installation and computational reproducibility. 

This pipeline currently only accepts paired-end reads as inputs. 

## Pipeline summary for metagenomic reads
1. Adapter trimming and quality control using fastp (0.22.0)
2. Removal of host (human) reads by mapping to a reference genome using bwa (0.7.17) 
3. Taxonomic classification of non-human reads using Kraken2 (2.1.2) and taxonomic abundance re-estimation using Bracken (2.6.1)

## Pipeline summary for metatranscriptomic reads
1. Adapter trimming and quality control using fastp (0.22.0)
2. Removal of host (human) reads using STAR (2.7.9a), a splice aware aligner.
3. Computational removal of prokaryotic and eukaryotic rRNAs using a k-mer based strategy with bbmap (38.93)
4. Optional sequence de-duplication using bbmap (38.93) clumpify.sh
5. Taxonomic classification of non-human reads using Kraken2 (2.1.2)
6. Functional annotation of reads by mapping to a microbial gene catalog of choice using bowtie2 (2.4.4)
7. Functional annotation of unmapped reads in the previous step, to a larger protein database e.g. [UHGP](http://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v1.0/uhgp_catalogue/) or [Uniref90](https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/) using diamond (2.0.12)
8. Pathway annotation of mappable reads using EggNOG mapper **WIP, unimplemented**

## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.04.0`) and add the nextflow executable to your $PATH

2. Install [`Docker`](https://docs.docker.com/engine/installation/)   

3. Clone the pipeline and refer to the help message
	```sh
	$ git clone https://github.com/Chiamh/meta-omics-nf
	
	$ nextflow run ./meta-omics-nf/main.nf --help
	```

	> Add a custom config file which contains the paths to various pre-installed databases. Refer to the test.config file in this repo for an example. 
	> Add a custom profile in the nextflow.config file, allowing you to specify the use of docker or singularity, and/or a task scheduler.  

4. Run the pipeline
	```sh
	$ nextflow run ./meta-omics-nf/main.nf -profile your_profile --rna_reads /path/to/metatranscriptomes --dna_reads /path/to/metagenomes --outdir /path/to/results
	```
	> You can specifiy multiple profiles separated by comma, e.g. -profile docker,sge.
	> The taxonomic classification, nucleotide alignment and translated search modules can be quite memory intensive depending on the databases used
	> Delete the work/ directory after running the pipeline to free up space taken up by intermediate files
	> There are modular workflows to reduce the size of intermediate files produced by the pipeline. See the help message for more details
	
## Contact

Minghao Chia: chia_minghao@gis.a-star.edu.sg, chiaminghao@gmail.com
