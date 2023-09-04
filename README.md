# Shotgun Metagenomics & Metatranscriptomics Nextflow Pipeline

## Introduction

Chiamh/meta-omics-nf is a bioinformatics pipeline that takes raw metagenomic **(MGX)** and/or metatranscriptomic **(MTX)** reads and annotates them at the taxonomic and functional levels.

The pipeline is built using [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) , a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. 
It uses Docker containers (also compatible with Singularity) for ease of installation and computational reproducibility. 

This pipeline currently only accepts paired-end reads as inputs. 

## Pipeline summary for metagenomic reads
1. Adapter trimming and quality control using fastp (0.22.0)
2. Optional removal of host (human) reads by mapping to a reference genome using bwa (0.7.17) 
3. Taxonomic classification of non-human reads using Kraken2 (2.1.2) and taxonomic abundance re-estimation using Bracken (2.6.1)
4. Optional removal of microbial spike-in sequences
5. Functional annotation of reads by mapping to a microbial gene catalog of choice using bowtie2 (2.4.4)
6. Functional annotation of unmapped reads in the previous step, to a larger protein database e.g. [UHGP](http://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v1.0/uhgp_catalogue/) or [Uniref90](https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/) using diamond (2.0.12)
7. Annotation of mappable reads into Clusters of Orthologous groups (COGs) using EggNOG mapper (2.1.6)

## Pipeline summary for metatranscriptomic reads
1. Adapter trimming and quality control using fastp (0.22.0)
2. Optional removal of host (human) reads using STAR (2.7.9a), a splice aware aligner.
3. Optional computational removal of prokaryotic and eukaryotic rRNAs using a k-mer based strategy with bbmap (38.93)
4. Optional sequence de-duplication using bbmap (38.93) clumpify.sh
5. Taxonomic classification of non-human reads using Kraken2 (2.1.2)
6. Functional annotation of reads by mapping to a microbial gene catalog of choice using bowtie2 (2.4.4)
7. Functional annotation of unmapped reads in the previous step, to a larger protein database e.g. [UHGP](http://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v1.0/uhgp_catalogue/) or [Uniref90](https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/) using diamond (2.0.12)
8. Annotation of mappable reads into Clusters of Orthologous groups (COGs) using EggNOG mapper (2.1.6)

<img src='/docs/metatrans_workflow.png' width='2000'>

## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.04.0`) and add the nextflow executable to your $PATH

2. Install [`Docker`](https://docs.docker.com/engine/installation/)   

3. Clone the pipeline and refer to the help message
	```sh
	$ git clone https://github.com/Chiamh/meta-omics-nf
	
	$ nextflow run ./meta-omics-nf/main.nf --help
	```
* Add a custom config file which contains the paths to various pre-installed databases. Refer to the test.config file in this repo for an example. 
* Add a custom profile in the nextflow.config file, allowing you to specify the use of docker or singularity, and/or a task scheduler.  

4. Make sure all helper scripts in meta-omics-nf/bin have execute permissions

	```sh
	$ chmod +x ./meta-omics-nf/bin/*
	```

5. Run the pipeline
	```sh
	$ nextflow run ./meta-omics-nf/main.nf -profile docker,your_profile --rna_reads /path/to/metatranscriptomes --dna_reads /path/to/metagenomes --outdir /path/to/results
	```
* You can specifiy multiple profiles separated by comma, e.g. -profile docker,sge.
* The taxonomic classification, nucleotide alignment, translated search and annotation modules can be quite memory intensive depending on the databases used
* Delete the work/ directory after running the pipeline to free up space taken up by intermediate files
* There are modular workflows (decontaminate and classify) to reduce the size of intermediate files produced by the pipeline. See the help message for more details.
* There is a concatenate workflow (-entry concatenate) to merge fastq.gz files across lanes for the same sample ID.
* You have the flexibility to turn off DNA spike in removal and/or the eggNOG annotation modules. See the help message for more details
	```sh
	$ nextflow run ./meta-omics-nf/main.nf -profile docker,your_profile -entry decontaminate --rna_reads /path/to/metatranscriptomes --dna_reads /path/to/metagenomes --outdir /path/to/results
	
	$ nextflow run ./meta-omics-nf/main.nf -profile docker,your_profile -entry classify --rna_reads /path/to/DECONTAMINATED_metatranscriptomes --dna_reads /path/to/DECONTAMINATED_metagenomes --outdir /path/to/results
	```
## Output files (bolded items are noteworthy for downstream analyses)

**Caveat: Read counts (paired-end) for some reports are not directly comparable with read counts (unpaired) for other reports. This is necessary because kraken2 is run in paired-end mode, whereas functional annotations take "unpaired" inputs.**

Why is it preferable to perform functional annotations using unpaired despite paired-end data? [Read this.](https://github.com/biobakery/humann#humann-30-and-paired-end-sequencing-data)

* decont/DNA (for metagenomes) or decont/RNA (for metatranscriptomes)
    * These folders contain decontaminated reads in fastq.gz format. Decontamination means adapter removal, host read removal, rRNA removal (for MTX) and de-duplication (for MTX). These steps are on by default, but are optional.

* kraken2_out/DNA
    * \*kraken2.out : Kraken2 [raw output](https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown#standard-kraken-output-format) with taxonomic classification for each read.
	* \*kraken2.tax : Kraken2 [report](https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown#sample-report-output-format)
	* \*kraken2_minimizer.tax : Kraken2 [report](https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown#distinct-minimizer-count-information) with distinct minimizer count information.
	* \*bracken.tax : Bracken report file with a simlar format to the Kraken2 report
	* **\*bracken.out** : Bracken output, tab separated with taxonomic classifications, read counts (paired-end) and abundance estimations.
	* **\*k2.s.tsv** : Simplified species level Kraken2 report, tab separated. Column 1: relative abundance; Column 2: paired read counts rooted at this taxon; Column 3: Number of minimizers in read data associated with this taxon; Column 4: An estimate of the number of distinct minimizers in read data associated with this taxon; Column 5: taxonomic classification at species level or unclassified. **Note: reads classified at genus level and above are not included in this report.**

* kraken2_out/RNA
    * \*kraken2.out : Kraken2 [raw output](https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown#standard-kraken-output-format) with taxonomic classification for each read.
	* \*kraken2.tax : Kraken2 [report](https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown#sample-report-output-format)
	* \*kraken2_minimizer.tax : Kraken2 [report](https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown#distinct-minimizer-count-information) with distinct minimizer count information.
	* **\*k2.s.tsv** : Simplified species level Kraken2 report, tab separated. Column 1: relative abundance; Column 2: paired read counts rooted at this taxon; Column 3: Number of minimizers in read data associated with this taxon; Column 4: An estimate of the number of distinct minimizers in read data associated with this taxon; Column 5: taxonomic classification at species level or unclassified. **Note: reads classified at genus level and above are not included in this report.**
	
* MGX_panalign_out/ (for metagenomes) or MTX_panalign_out/ (for metatranscriptomes)
    * \*bt2_pangenome_aligned.bam : BAM file after alignment of reads (unpaired) to pangene catalog.
	* \*bt2_pangenome_aligned_filtered_cov.tsv : Tab separated file containing read coverage across pangenes. **Only pangenes with >= 50% coverage are reported here.**
	* \*bt2_pangenome_unaligned.fastq.gz : All reads that did not align to the pangene catalog.

* MGX_dmnd_out/ (for metagenomes) or MTX_panalign_out/ (for metatranscriptomes)
    * \*uniref90_aligned.out : Tab separated file containing alignment information to Uniref90 database. Columns are in [BLAST format](https://www.metagenomics.wiki/tools/blast/blastn-output-format-6) in this order: qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen evalue bitscore salltitles
	* \*uniref90_unaligned.fa : All reads that did not align to the Uniref90 database (and also by consequence, the pangene catalog)
	
* MGX_annotations/ (for metagenomes) or MTX_annotations/ (for metatranscriptomes)
    * \*decont.emapper.annotations : Tab separated output from eggnog mapper 2.1.6.\
**Note: Apart from reads which did not align to pangene catalogs or Uniref90 databases, the following annotations are over features with >= 50% read coverage.**
	* **\*all_aligned_taxonomy.tsv** : Tab separated output. Column 1: Pangene ID; Column 2: Read ID; Column 3: Kraken2 taxonomic classification; Column 4: Uniref90 ID.
	* **\*all_aligned_taxonomy_summary.tsv** : Tab separated output summarizing read counts per kraken2 taxon per gene feature. Column 1: Number of unpaired reads; Column 2: Pangene ID; Column 3: Kraken2 classification; Column 4: Uniref90 ID.
	* **\*unaligned_taxonomy.tsv** : Tab separated output. Column 1: Read ID; Column 2: Kraken2 taxonomic classification. This details the taxonomic classifications for reads which did not align to pangene catalogs or Uniref90 databases
	* **\*panalign_annot.tsv** : Tab separated output containing functional annotations (COGs, Uniref90 and Gene Onthologies) for pangenes with >= 50% read coverage. Columns are in order: pangene length, percent_cov, unpaired_read_count, pangene_desc, uniref90_ID, uniref90_desc, uniref90_GO, emapper_max_annot_OG, emapper_OG, [emapper_max_annot_lvl](https://github.com/eggnogdb/eggnog-mapper/wiki/eggNOG-mapper-v2.1.5-to-v2.1.12#user-content-Annotations_file), [eggnog_cat](http://clovr.org/docs/clusters-of-orthologous-groups-cogs/), eggnog_GO, eggnog_desc
	* **\*transl-search_annot.tsv** : Tab separated output containing functional annotations (COGs, Uniref90 and Gene Onthologies) for Uniref90 clusters (with no pangene assignment) with >= 50% read coverage. Columns are in order: uniref90_ID, percent_cov, AA_length, unpaired_read_count, uniref90_desc, emapper_max_annot_OG, emapper_OG, emapper_max_annot_lvl, eggnog_cat, eggnog_GO, eggnog_desc, uniref90_GO

## Updates

* There is now a "concatenate" subworkflow to merge fastq.gz files across different lanes for the same sample ID.
* Fixed code to not have the same variable name appear more than once in the input block for any process. 
* Add distinct minimizer information to kraken2 reports for false positive filtering. (Mar 2023)
* Add taxonomy summary reports (Mar 2023)

	
## Contact

Minghao Chia: chia_minghao@gis.a-star.edu.sg, chiaminghao@gmail.com
