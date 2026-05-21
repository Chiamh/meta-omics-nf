

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


if (!params.human_pangenome_path && !params.decont_off && params.process_dna){
    helpMessage()
    log.info"""
    [Error] --human_pangenome_path is required for removal of host (human) reads from metagenomic sequences (decontamination step)
    """.stripIndent()
    exit 0
}

if (!params.star_index && !params.decont_off){
    helpMessage()
    log.info"""
    [Error] --star_index is required for removal of host (human) reads (decontamination step)
    """.stripIndent()
    exit 0
}

if (!params.ribokmers && !params.decont_off && params.process_rna && params.remove_rRNA){
    helpMessage()
    log.info"""
    [Error] --ribokmers is required for removal of rRNA reads (decontamination steps)
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
        def read1_files = files("${params.rna_reads}/**/*${row.id}*{_R1,_1}.{fastq,fq}.gz") +
                          files("${params.rna_reads}/*${row.id}*{_R1,_1}.{fastq,fq}.gz")
        
        def read2_files = files("${params.rna_reads}/**/*${row.id}*{_R2,_2}.{fastq,fq}.gz") + 
                          files("${params.rna_reads}/*${row.id}*{_R2,_2}.{fastq,fq}.gz")
        
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
        def read1_files = files("${params.dna_reads}/**/*${row.id}*{_R1,_1}.{fastq,fq}.gz") +
                          files("${params.dna_reads}/*${row.id}*{_R1,_1}.{fastq,fq}.gz")
        
        def read2_files = files("${params.dna_reads}/**/*${row.id}*{_R2,_2}.{fastq,fq}.gz") + 
                          files("${params.dna_reads}/*${row.id}*{_R2,_2}.{fastq,fq}.gz")
        
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

include { FASTP_UMI } from '../modules/fastp_umi.nf'
include { STAR } from '../modules/star.nf'
include { RIBOFILTER } from '../modules/rRNAfilter.nf'
include { DECONT_DNA } from '../modules/decont_dna.nf'
include { DECONT_DNA_PANALIGN } from '../modules/decont_dna_panalign.nf'

/*
========================================================================================
    Named workflow
========================================================================================
*/

workflow DECONT {

if (params.process_rna){
     if ( params.process_rna ){
        // 1. Read QC and add umi to headers using fastp
        FASTP_UMI( ch_rna_input )
        ch_rna_input = FASTP_UMI.out.reads

        // 2. remove human host RNA and 3. remove rRNAs
        if ( params.decont_off ){
            if ( params.remove_rRNA ){
                RIBOFILTER(params.ribokmers, ch_rna_input)
                ch_rna_decont = RIBOFILTER.out.reads
            } else {
                ch_rna_decont = ch_rna_input
            }
        } else {
            STAR(params.star_index, ch_rna_input)
            if ( params.remove_rRNA ){
                RIBOFILTER(params.ribokmers, STAR.out.microbereads)
                ch_rna_decont = RIBOFILTER.out.reads
            } else {
                ch_rna_decont = STAR.out.microbereads
            }
        }
	}
}

if (params.process_dna){
     if ( params.decont_off ) {
            ch_dna_decont = ch_dna_input
        } else {
            DECONT_DNA(params.star_index, ch_dna_input)
            DECONT_DNA_PANALIGN(params.human_pangenome_path, DECONT_DNA.out.reads)
            ch_dna_decont = DECONT_DNA_PANALIGN.out.reads
        }
}
}






