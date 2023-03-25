#!/usr/bin/env nextflow

/*
========================================================================================
    Shotgun metagenomics and metatranscriptomics 
========================================================================================
    Github : https://github.com/Chiamh/meta-omics-nf
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl=2

/*
========================================================================================
    Help messages and warnings
========================================================================================
*/

//https://www.baeldung.com/groovy-def-keyword
//https://www.nextflow.io/blog/2020/cli-docs-release.html

def helpMessage() {
  // adapted from nf-core
    log.info"""
    
    Usage for main workflow:
    The typical command for running the pipeline is as follows:
      nextflow run main.nf
      --rna_reads FOLDER_FOR_RNA_READS
      --dna_reads FOLDER_FOR_DNA_READS
      --bwaidx_path FOLDER_FOR_HUMAN_GENOME_AND_BWA_INDEX
      --bwaidx NAME_OF_BWA_INDEX
      --star_index FOLDER_FOR_STAR_INDEX_FOR_HUMAN_GENOME
      --ribokmers FOLDER_FOR_BBMAP_RIBOKMERS
      --kraken2db FOLDER_FOR_KRAKEN2_AND_BRACKEN_DB
      --pangenome_path FOLDER_FOR_PANGENOME_AND_BOWTIE2_INDEX
	  --pangenome	NAME_OF_PANGENOME_BOWTIE2_INDEX
      --dmnddb PATH_TO_DIAMOND2_DB
    
	NOTE: A more user-friendly approach is to specify these parameters in a *.config file under a custom profile 
	
    IMPT: Set either the --process_rna or --process_dna arguments to false if no RNA or DNA reads are provided, respsectively. 
    
    The main workflow can take up a lot of disk space with intermediate fastq files. 
    
    If this is a problem, the workflow can be run as two separate modules.
    The decontaminate module removes the intermediate files and is not compatible with nextflow -resume
    The classify module is still compatible with nextflow -resume because the smaller intemediate files are kept in the nextflow work/ directory
    The classify module assumes gzipped compressed reads as direct inputs to Kraken2
    
    Usage for alternative workflow:
    nextflow run main.nf -entry decontaminate [args]...
    nextflow run main.nf -entry classify [args]...
    
    Input and database arguments are null by default.
    Rather than manually specifying the paths to so many databases, it is best to create a custom nextflow config file.
     
    Input arguments:
      --rna_reads                   Path to a folder containing all input metatranscriptomic reads (this will be recursively searched for *fastq.gz/*fq.gz/*fq/*fastq files)
      --dna_reads                   Path to a folder containing all input metagenomic reads (this will be recursively searched for *fastq.gz/*fq.gz/*fq/*fastq files)
    Database arguments:
      --bwaidx_path                 Path to the folder with host (human) reference genome and bwa index
      --bwaidx			    Name of the bwa index e.g. hg38.fa
      --star_index                  Path to the directory containing the index for the human genome for STAR aligner
      --ribokmers                   Path to the eukaryotic and prokaryotic ribokmer database for computational rRNA removal using BBmap
      --kraken2db                   Path to the Kraken2 and Bracken databases
      --pangenome_path              Path to the folder with bowtie2 index for custom-built microbial pangenome/gene catalog
      --pangenome                   Name of the bowtie2 index for the pangenome/gene catalog e.g. IHSMGC
      --dmnddb                      Path to a custom-built Diamond 2 database (e.g. *.dmnd)
      --eggnog_db                   Path to folder containing the eggnog database
      --eggnog_OG_annots            Path to a pre-built e5.og_annotations.tsv file, downloaded from http://eggnog5.embl.de/download/eggnog_5.0, sorted by EGGNOG ID
      --uniref90_fasta              Path to fasta file containing amino acid sequences from Uniref90
      --uniref90_GO                 Path to two column .tsv file derived from https://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/idmapping/idmapping_selected.tab.gz
      --pangenome_annots            Path to pre-computed eggnog annotations for pangenome
      --spike_in_path		    Path to file denoting genera/species to remove from metagenomes before functional profiling, because they are spike-ins
    Bracken options:
      --readlength                  Length of Bracken k-mers to use [default: 150]
    Workflow options:
      -entry                        Can be one of [decontaminate, classify]. For disk space saving workflows. Note SINGLE dash.
      --process_rna                 Turns on steps to process metatranscriptomes [Default: true]. If true, --rna_reads is a mandatory argument
      --process_dna                 Turns on steps to process metagenomes [Default: true]. If true, --dna_reads is a mandatory argument
      --decont_off                  Skip trimming, QC and decontamination steps [Default: false]
      --profilers_off               Skip Kraken2 and Bracken steps [Default: false]
      --panalign_off                Skip pangenome alignment with bowtie 2. Will also skip translated search with Diamond [Default: false]
      --diamond_off                 Skip translated search with Diamond [Default: false]
      --rm_spikes                   Removes spike in sequences from metagenomes [Default: true]
      --annotate_off                Skip functional annotation using Eggnog and Uniref90 [Default: false]
    Output arguments:
      --outdir                      The output directory where the results will be saved [Default: ./pipeline_results]
      --tracedir                    The directory where nextflow logs will be saved [Default: ./pipeline_results/pipeline_info]
    AWSBatch arguments:
      --awsregion                   The AWS Region for your AWS Batch job to run on [Default: false]
      --awsqueue                    The AWS queue for your AWS Batch job to run on [Default: false]
    Others:
      --help		            Display this help message
    """
}

if (params.help){
    helpMessage()
    exit 0
}


/*
========================================================================================
    Include modules
========================================================================================
*/

include { FULL } from './workflows/full_workflow.nf'
include { DECONT } from './workflows/decont.nf'
include { PROFILE } from './workflows/classify.nf'
include { CONCATENATE } from './workflows/concatenate.nf'

/*
========================================================================================
    Main workflow (default)
========================================================================================
*/

// this main workflow will generate all intermediate files and can be resumed with nextflow
workflow {
    
    FULL ()
     
}

//modular named workflows to reduce intermediate file size. 
//warning: the decontaminate named workflow is not compatible with nextflow -resume because intermediate files are deleted to save space
workflow decontaminate {
    
    DECONT ()
}

// This classify workflow can be resumed with nextflow
// Typical use is to specify the path to already decontaminated reads with --rna_reads and/or --dna_reads
// Assumes gzipped compressed reads
workflow classify {
    
    PROFILE ()
}
// Use the concatenate workflow to join fastq files across different lanes by libid. Specify path to raw reads with --rna_reads and/or --dna_reads
workflow concatenate {

    CONCATENATE ()
}
