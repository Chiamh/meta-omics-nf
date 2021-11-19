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
    Bracken options:
      --readlength                  Length of Bracken k-mers to use [default: 150]
    Workflow options:
      -entry                        Can be one of [decontaminate, classify]. For disk space saving workflows. Note SINGLE dash.
      --process_rna                 Turns on steps to process metatranscriptomes [Default: true]. If true, --rna_reads is a mandatory argument
      --process_dna                 Turns on steps to process metagenomes [Default: true]. If true, --dna_reads is a mandatory argument
      --decont_off                  Skip trimming, QC and decontamination steps [Default: false]
      --dedupe                      Run sequence deduplication with clumpify.sh [Default: true]
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
      --help		            Display this help message
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

if (!params.kraken2db && !params.profilers_off){
    helpMessage()
    log.info"""
    [Error] --kraken2db is required for taxonomic classification
    """.stripIndent()
    exit 0
}

if (!params.pangenome_path && !params.panalign_off){
    helpMessage()
    log.info"""
    [Error] --pangenome_path is required for mapping of metatranscriptomes to gene catalog
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
    Define channels for read pairs
========================================================================================
*/

//params.rna_reads = "$baseDir/data/raw_fastq/rna/*{1,2}.{fq.gz,fastq.gz}"
//params.dna_reads = "$baseDir/data/raw_fastq/dna/*{1,2}.{fq.gz,fastq.gz}"

if (params.process_rna){
    Channel.fromFilePairs( [params.rna_reads + '/**{R,.,_}{1,2}*{fastq,fastq.gz,fq,fq.gz}'], checkExists:true ).set{ ch_rna_input }
}

if (params.process_dna){
    Channel.fromFilePairs( [params.dna_reads + '/**{R,.,_}{1,2}*{fastq,fastq.gz,fq,fq.gz}'], checkExists:true ).set{ ch_dna_input }
}

/*
========================================================================================
    Include modules
========================================================================================
*/

include { KRAKEN2_RNA } from '../modules/kraken_rna.nf'
include { KRAKEN2_DNA } from '../modules/kraken_dna.nf'
include { BRACKEN } from '../modules/bracken.nf'
include { PANALIGN } from '../modules/panalign.nf'
include { DMND } from '../modules/dmnd.nf'

/*
========================================================================================
    Named workflow
========================================================================================
*/

workflow PROFILE {
if (params.process_rna){
   KRAKEN2_RNA(params.kraken2db, ch_rna_input)
   PANALIGN(params.pangenome_path, ch_rna_input)
   DMND(params.dmnddb, PANALIGN.out.unaligned)
}
if (params.process_dna){   
   KRAKEN2_DNA(params.kraken2db, ch_dna_input)
   BRACKEN(params.kraken2db, params.readlength, KRAKEN2_DNA.out.k2tax)
}
}

