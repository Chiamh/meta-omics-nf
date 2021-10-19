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
      --hg_fasta FOLDER_FOR_HUMAN_GENOME_AND_BWA_INDEX
      --star_index FOLDER_FOR_STAR_INDEX_FOR_HUMAN_GENOME
      --ribokmers FOLDER_FOR_BBMAP_RIBOKMERS
      --kraken2db FOLDER_FOR_KRAKEN2_AND_BRACKEN_DB
      --pangenome FOLDER_FOR_PANGENOME_AND_BOWTIE2_INDEX
      --dmnddb FOLDER_FOR_DIAMOND2_DB
    
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
      --hg_fasta                    Path to the host (human) reference genome with bwa index in the same folder
      --star_index                  Path to the directory containing the index for the human genome for STAR aligner
      --ribokmers                   Path to the eukaryotic and prokaryotic ribokmer database for computational rRNA removal using BBmap
      --kraken2db                   Path to the Kraken2 and Bracken databases
      --pangenome                   Path to a custom-built microbial pangenome with bowtie2 index
      --dmnddb                      Path to a custom-build Diamond 2 database (usually Uniref90)  
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
    Output arguments:
      --outdir                      The output directory where the results will be saved [Default: ./pipeline_results]
      --tracedir                    The directory where nextflow logs will be saved [Default: ./pipeline_results/pipeline_info]
    AWSBatch arguments:
      --awsregion                   The AWS Region for your AWS Batch job to run on [Default: false]
      --awsqueue                    The AWS queue for your AWS Batch job to run on [Default: false]
    Others:
      --help 			    Display this help message
    """
}

if (params.help){
    helpMessage()
    exit 0
}


if (!params.rna_reads && !params.dna_reads){
    helpMessage()
    log.info"""
    [Error] The path to at least one input folder for sequences is required
    """.stripIndent()
    exit 0
}

if (!params.rna_reads && params.process_rna){
    helpMessage()
    log.info"""
    [Error] The path to input RNA sequences is required because --process_rna is true
    """.stripIndent()
    exit 0
}

if (!params.dna_reads && params.process_dna){
    helpMessage()
    log.info"""
    [Error] The path to input DNA sequences is required because --process_dna is true
    """.stripIndent()
    exit 0
}


if (!params.hg_fasta && !params.decont_off && params.process_dna){
    helpMessage()
    log.info"""
    [Error] --hg_fasta is required for removal of host (human) reads from metagenomic sequences (decontamination step)
    """.stripIndent()
    exit 0
}

if (!params.star_index && !params.decont_off && params.process_rna){
    helpMessage()
    log.info"""
    [Error] --star_index is required for removal of host (human) reads from metatranscriptomic sequences (decontamination step)
    """.stripIndent()
    exit 0
}

if (!params.ribokmers && !params.decont_off && params.process_rna){
    helpMessage()
    log.info"""
    [Error] --ribokmers is required for removal of rRNA reads (decontamination steps)
    """.stripIndent()
    exit 0
}

if (!params.kraken2db && !params.profilers_off){
    helpMessage()
    log.info"""
    [Error] --kraken2db is required for taxonomic classification
    """.stripIndent()
    exit 0
}

if (!params.pangenome && !params.panalign_off){
    helpMessage()
    log.info"""
    [Error] --pangenome is required for mapping of metatranscriptomes to gene catalog
    """.stripIndent()
    exit 0
}

if (!params.dmnddb && !params.diamond_off){
    helpMessage()
    log.info"""
    [Error] --dmnddb is required for translated search
    """.stripIndent()
    exit 0
}

/*
========================================================================================
    Include modules
========================================================================================
*/

include { FASTP } from '../modules/fastp.nf'
include { RIBOFILTER } from '../modules/rRNAfilter.nf'
include { DEDUP } from '../modules/dedup.nf'
include { KRAKEN2_RNA } from '../modules/kraken_rna.nf'
include { KRAKEN2_DNA } from '../modules/kraken_dna.nf'
include { BRACKEN } from '../modules/bracken.nf'
include { PANALIGN } from '../modules/panalign.nf'
include { DMND } from '../modules/dmnd.nf'
include { DECONT_DNA } from '../modules/decont_dna.nf'

/*
========================================================================================
    Workflow
========================================================================================
*/

//https://bioinformatics.stackexchange.com/questions/15739/use-conditional-in-workflow-in-nextflow-dsl2
//https://github.com/nextflow-io/patterns/blob/master/docs/conditional-process.adoc

// this main workflow will generate all intermediate files and can be resumed with nextflow
workflow FULL {
    if ( params.process_rna ){
        FASTP(params.star_index, ch_rna_input)
        RIBOFILTER(params.ribokmers, FASTP.out.microbereads)
        DEDUP(RIBOFILTER.out.reads)
        
        if ( params.decont_off ) {
            ch_rna_decont = ch_rna_input
        } else if ( !params.decont_off && params.dedupe ){
            ch_rna_decont = DEDUP.out.reads
            } else if ( !params.decont_off && !params.dedupe ){
            ch_rna_decont = RIBOFILTER.out.reads
    }
        KRAKEN2_RNA(params.kraken2db, ch_rna_decont)
        PANALIGN(params.pangenome, ch_rna_decont)
        DMND(params.dmnddb, PANALIGN.out.unaligned)
    }
    
    if ( params.process_dna ){
        DECONT_DNA(params.hg_fasta, ch_dna_input)
        
        if ( params.decont_off ) {
            ch_dna_decont = ch_dna_input
        } else if ( !params.decont_off ) {
            ch_dna_decont = DECONT_DNA.out.reads
        }
        
        KRAKEN2_DNA(params.kraken2db, ch_dna_decont)
        BRACKEN(params.kraken2db, params.readlength, KRAKEN2_DNA.out.k2tax)
    }  
}






