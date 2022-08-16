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
      --eggnog_db		    Path to folder containing the eggnog database
      --eggnog_OG_annots	    Path to a pre-built e5.og_annotations.tsv file, downloaded from http://eggnog5.embl.de/download/eggnog_5.0, sorted by EGGNOG ID
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
      --rm_spikes		    Removes spike in sequences from metagenomes [Default: true]
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


if (!params.bwaidx_path && !params.decont_off && params.process_dna){
    helpMessage()
    log.info"""
    [Error] --bwaidx_path is required for removal of host (human) reads from metagenomic sequences (decontamination step)
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

if (!params.pangenome_path && !params.panalign_off){
    helpMessage()
    log.info"""
    [Error] --pangenome_path is required for mapping of MTX/MGX reads to gene catalog
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

if (!params.eggnog_db && !params.annotate_off){
    helpMessage()
    log.info"""
    [Error] --eggnog_db is required for functional annotation. This is the path to the folder containing the eggnog database.
    """.stripIndent()
    exit 0
}

if (!params.eggnog_OG_annots && !params.annotate_off){
    helpMessage()
    log.info"""
    [Error] --eggnog_OG_annots is required for functional annotation. This is the path to a pre-built e5.og_annotations.tsv file, downloaded from http://eggnog5.embl.de/download/eggnog_5.0, sorted by EGGNOG ID
    """.stripIndent()
    exit 0
}

if (!params.uniref90_fasta && !params.annotate_off){
    helpMessage()
    log.info"""
    [Error] --uniref90_fasta is required for functional annotation. This is a fasta file containing amino acid sequences from Uniref90
    """.stripIndent()
    exit 0
}

if (!params.uniref90_GO && !params.annotate_off){
    helpMessage()
    log.info"""
    [Error] --uniref90_GO is required for functional annotation. This is a two column .tsv file derived from https://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/idmapping/idmapping_selected.tab.gz.
    """.stripIndent()
    exit 0
}

if (!params.pangenome_annots && !params.annotate_off){
    helpMessage()
    log.info"""
    [Error] --pangenome_annots is required for functional annotation. Check documentation for the format of this file
    """.stripIndent()
    exit 0
}

if (!params.spike_in_path && params.rm_spikes){
   helpMessage()
   log.info"""
   [Error] --spike_in_path is required to remove spike ins from MGX data. This is the path to a file containing the spike in genera/species to remove
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
    Channel.fromFilePairs( [params.rna_reads + '/**{R,.,_}{1,2}*{fastq,fastq.gz,fq,fq.gz}'], checkIfExists:true ).set{ ch_rna_input }
}

if (params.process_dna){
    Channel.fromFilePairs( [params.dna_reads + '/**{R,.,_}{1,2}*{fastq,fastq.gz,fq,fq.gz}'], checkIfExists:true ).set{ ch_dna_input }
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
include { PANALIGN_RNA } from '../modules/panalign_rna.nf'
include { PANALIGN_DNA } from '../modules/panalign_dna.nf'
include { PANALIGN_DNA_SPIKES } from '../modules/panalign_dna_spikes.nf'
include { DMND_RNA } from '../modules/dmnd_rna.nf'
include { DMND_DNA } from '../modules/dmnd_dna.nf'
include { DECONT_DNA } from '../modules/decont_dna.nf'
include { ANNOT_DMND_RNA } from '../modules/annot_dmnd_rna.nf'
include { ANNOT_DMND_DNA } from '../modules/annot_dmnd_dna.nf'
include { ANNOT_PAN_RNA } from '../modules/annot_pan_rna.nf'
include { ANNOT_PAN_DNA } from '../modules/annot_pan_dna.nf'
include { TRF_TAXA_DNA } from '../modules/transfer_taxa_dna.nf'
include { TRF_TAXA_RNA } from '../modules/transfer_taxa_rna.nf'

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
	
	if ( params.dedupe ){
        DEDUP(RIBOFILTER.out.reads)
	}
        
        if ( params.decont_off ) {
            ch_rna_decont = ch_rna_input
        } else if ( !params.decont_off && params.dedupe ){
            ch_rna_decont = DEDUP.out.reads
            } else if ( !params.decont_off && !params.dedupe ){
            ch_rna_decont = RIBOFILTER.out.reads
    }
        KRAKEN2_RNA(params.kraken2db, ch_rna_decont)
        PANALIGN_RNA(params.pangenome_path, ch_rna_decont)
        DMND_RNA(params.dmnddb, PANALIGN_RNA.out.unaligned)
	ANNOT_DMND_RNA(params.uniref90_fasta, params.eggnog_OG_annots, params.eggnog_db, params.uniref90_GO, DMND_RNA.out.aligned)
	ANNOT_PAN_RNA(params.pangenome_annots, PANALIGN_RNA.out.coverage)

	ch_trf_taxa_rna_in=KRAKEN2_RNA.out.k2out.join(PANALIGN_RNA.out.aligned).join(DMND_RNA.out.aligned).join(DMND_RNA.out.unaligned)
	TRF_TAXA_RNA(params.pangenome_annots, ch_trf_taxa_rna_in)
    }
    
    if ( params.process_dna ){
        DECONT_DNA(params.bwaidx_path, ch_dna_input)
        
        if ( params.decont_off ) {
            ch_dna_decont = ch_dna_input
        } else if ( !params.decont_off ) {
            ch_dna_decont = DECONT_DNA.out.reads
        }
        
        KRAKEN2_DNA(params.kraken2db, ch_dna_decont)
        BRACKEN(params.kraken2db, params.readlength, KRAKEN2_DNA.out.k2tax)

	if ( params.rm_spikes ){

	ch_panalign_dna_in=ch_dna_decont.join(KRAKEN2_DNA.out.k2out)

	PANALIGN_DNA_SPIKES(params.pangenome_path, params.spike_in_path, ch_panalign_dna_in)
	DMND_DNA(params.dmnddb, PANALIGN_DNA_SPIKES.out.unaligned)
	ANNOT_DMND_DNA(params.uniref90_fasta, params.eggnog_OG_annots, params.eggnog_db, params.uniref90_GO, DMND_DNA.out.aligned)
	ANNOT_PAN_DNA(params.pangenome_annots, PANALIGN_DNA_SPIKES.out.coverage)

	ch_trf_taxa_dna_in=KRAKEN2_DNA.out.k2out.join(PANALIGN_DNA_SPIKES.out.aligned).join(DMND_DNA.out.aligned).join(DMND_DNA.out.unaligned)

	TRF_TAXA_DNA(params.pangenome_annots, ch_trf_taxa_dna_in)
	} else if ( !params.rm_spikes ){
	PANALIGN_DNA(params.pangenome_path, ch_dna_decont)
	DMND_DNA(params.dmnddb, PANALIGN_DNA.out.unaligned)
	ANNOT_DMND_DNA(params.uniref90_fasta, params.eggnog_OG_annots, params.eggnog_db, params.uniref90_GO, DMND_DNA.out.aligned)
	ANNOT_PAN_DNA(params.pangenome_annots, PANALIGN_DNA.out.coverage)

	ch_trf_taxa_dna_in=KRAKEN2_DNA.out.k2out.join(PANALIGN_DNA.out.aligned).join(DMND_DNA.out.aligned).join(DMND_DNA.out.unaligned)

	TRF_TAXA_DNA(params.pangenome_annots, ch_trf_taxa_dna_in)
	}	
    }  
}
