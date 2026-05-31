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
    [Error] --pangenome_path is required for mapping of reads to pangenome or gene catalog
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

if (params.process_rna && params.rna_list){
        Channel
        .fromPath(params.rna_list)
        .splitCsv(header: true)
        .map { row ->
        // Recursively find files matching the sample name pattern
        def read1_files = files("${params.rna_reads}/**/*${row.id}*_R1_*.{fastq,fq}.gz") +
                          files("${params.rna_reads}/*${row.id}*_R1_*.{fastq,fq}.gz") +
                          files("${params.rna_reads}/**/*${row.id}*_1.{fastq,fq}.gz") +
                          files("${params.rna_reads}/*${row.id}*_1.{fastq,fq}.gz")
						  
        def read2_files = files("${params.rna_reads}/**/*${row.id}*_R2_*.{fastq,fq}.gz") +
                          files("${params.rna_reads}/*${row.id}*_R2_*.{fastq,fq}.gz") +
                          files("${params.rna_reads}/**/*${row.id}*_2.{fastq,fq}.gz") +
                          files("${params.rna_reads}/*${row.id}*_2.{fastq,fq}.gz")
        
        read1_files = read1_files.unique()
        read2_files = read2_files.unique()
        
        if (read1_files.size() == 0) error "No R1 found for RNA: ${row.id}"
        if (read2_files.size() == 0) error "No R2 found for RNA: ${row.id}"
        

        tuple(row.id, tuple(read1_files[0], read2_files[0]))
        }
        .set{ ch_rna_input }
} else if (params.process_rna && !params.rna_list){
        Channel.fromFilePairs( [params.rna_reads + '/**{R,.,_}{1,2}*{fastq,fastq.gz,fq,fq.gz}'], checkIfExists:true ).set{ ch_rna_input }
}

if (params.process_dna && params.dna_list){
        Channel
        .fromPath(params.dna_list)
        .splitCsv(header: true)
        .map { row ->
        // Recursively find files matching the sample name pattern
        def read1_files = files("${params.dna_reads}/**/*${row.id}*_R1_*.{fastq,fq}.gz") +
                          files("${params.dna_reads}/*${row.id}*_R1_*.{fastq,fq}.gz") +
                          files("${params.dna_reads}/**/*${row.id}*_1.{fastq,fq}.gz") +
                          files("${params.dna_reads}/*${row.id}*_1.{fastq,fq}.gz")
						  
        def read2_files = files("${params.dna_reads}/**/*${row.id}*_R2_*.{fastq,fq}.gz") +
                          files("${params.dna_reads}/*${row.id}*_R2_*.{fastq,fq}.gz") +
                          files("${params.dna_reads}/**/*${row.id}*_2.{fastq,fq}.gz") +
                          files("${params.dna_reads}/*${row.id}*_2.{fastq,fq}.gz")
        
        read1_files = read1_files.unique()
        read2_files = read2_files.unique()
        
        if (read1_files.size() == 0) error "No R1 found for DNA: ${row.id}"
        if (read2_files.size() == 0) error "No R2 found for DNA: ${row.id}"
        

        tuple(row.id, tuple(read1_files[0], read2_files[0]))
        }
        .set{ ch_dna_input }
} else if (params.process_dna && !params.dna_list){
        Channel.fromFilePairs( [params.dna_reads + '/**{R,.,_}{1,2}*{fastq,fastq.gz,fq,fq.gz}'], checkIfExists:true ).set{ ch_dna_input }
}

/*
========================================================================================
    Include modules
========================================================================================
*/

include { KRAKEN2_RNA } from '../modules/kraken_rna.nf'
include { KRAKEN2_DNA } from '../modules/kraken_dna.nf'
include { BRACKEN } from '../modules/bracken.nf'
include { PANALIGN_RNA } from '../modules/panalign_rna.nf'
include { PANALIGN_DNA } from '../modules/panalign_dna.nf'
include { PANALIGN_DNA_SPIKES } from '../modules/panalign_dna_spikes.nf'
include { DMND_RNA } from '../modules/dmnd_rna.nf'
include { DMND_DNA } from '../modules/dmnd_dna.nf'
include { ANNOT_DMND_RNA } from '../modules/annot_dmnd_rna.nf'
include { ANNOT_DMND_DNA } from '../modules/annot_dmnd_dna.nf'
include { ANNOT_PAN_RNA } from '../modules/annot_pan_rna.nf'
include { ANNOT_PAN_DNA } from '../modules/annot_pan_dna.nf'
include { TRF_TAXA_DNA } from '../modules/transfer_taxa_dna.nf'
include { TRF_TAXA_RNA } from '../modules/transfer_taxa_rna.nf'

/*
========================================================================================
    Named workflow
========================================================================================
*/

workflow PROFILE {
if (params.process_rna){
   KRAKEN2_RNA(params.kraken2db, ch_rna_input)
   PANALIGN_RNA(params.pangenome_path, ch_rna_input)
   DMND_RNA(params.dmnddb, PANALIGN_RNA.out.unaligned)
   ANNOT_DMND_RNA(params.uniref90_fasta, params.eggnog_OG_annots, params.eggnog_db, params.uniref90_GO, DMND_RNA.out.aligned)
   ANNOT_PAN_RNA(params.pangenome_annots, PANALIGN_RNA.out.coverage)

   ch_trf_taxa_rna_in=KRAKEN2_RNA.out.k2out.join(PANALIGN_RNA.out.aligned).join(DMND_RNA.out.aligned).join(DMND_RNA.out.unaligned)
   TRF_TAXA_RNA(params.pangenome_annots, ch_trf_taxa_rna_in)
}
if (params.process_dna){
   KRAKEN2_DNA(params.kraken2db, ch_dna_input)
   BRACKEN(params.kraken2db, params.readlength, KRAKEN2_DNA.out.k2tax)
   if ( params.rm_spikes ){

   ch_panalign_dna_in=ch_dna_input.join(KRAKEN2_DNA.out.k2out)
   PANALIGN_DNA_SPIKES(params.pangenome_path, params.spike_in_path, ch_panalign_dna_in)
   DMND_DNA(params.dmnddb, PANALIGN_DNA_SPIKES.out.unaligned)
   ANNOT_DMND_DNA(params.uniref90_fasta, params.eggnog_OG_annots, params.eggnog_db, params.uniref90_GO, DMND_DNA.out.aligned)
   ANNOT_PAN_DNA(params.pangenome_annots, PANALIGN_DNA_SPIKES.out.coverage)

   ch_trf_taxa_dna_in=KRAKEN2_DNA.out.k2out.join(PANALIGN_DNA_SPIKES.out.aligned).join(DMND_DNA.out.aligned).join(DMND_DNA.out.unaligned)
   TRF_TAXA_DNA(params.pangenome_annots, ch_trf_taxa_dna_in)
   } else if ( !params.rm_spikes ){
   PANALIGN_DNA(params.pangenome_path, ch_dna_input)
   DMND_DNA(params.dmnddb, PANALIGN_DNA.out.unaligned)
   ANNOT_DMND_DNA(params.uniref90_fasta, params.eggnog_OG_annots, params.eggnog_db, params.uniref90_GO, DMND_DNA.out.aligned)
   ANNOT_PAN_DNA(params.pangenome_annots, PANALIGN_DNA.out.coverage)

   ch_trf_taxa_dna_in=KRAKEN2_DNA.out.k2out.join(PANALIGN_DNA.out.aligned).join(DMND_DNA.out.aligned).join(DMND_DNA.out.unaligned)
   TRF_TAXA_DNA(params.pangenome_annots, ch_trf_taxa_dna_in)
   }
}
}
