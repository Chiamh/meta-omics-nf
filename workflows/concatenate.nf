// WORKFLOW file for concatenating R1 and R2 fastq.gz files that are spread across different lanes, using a common library ID
// Assumes that the R1 and R2 files are in the same subdirectory, named after a LIBID like MHS584. e.g. $INPUTDIR/MHS584/*
// params.rna_reads and params.dna_reads are paths to folders containing these subdirectories
// Yes you could create a glob pattern parameter to make this more generalizable
/*
========================================================================================
    Define channels for read pairs
========================================================================================
*/

if (params.process_rna){
    rna_reads = Channel.fromPath( [params.rna_reads + '/**{R,.,_}{1,2}*{fastq,fastq.gz,fq,fq.gz}'], checkIfExists:true ).map{
								file -> tuple( file.getParent().getName(), file )}.groupTuple(sort:true)
}

if (params.process_dna){
    dna_reads = Channel.fromPath( [params.dna_reads + '/**{R,.,_}{1,2}*{fastq,fastq.gz,fq,fq.gz}'], checkIfExists:true ).map{
								file -> tuple( file.getParent().getName(), file )}.groupTuple(sort:true)	
}

/*
========================================================================================
    Include modules
========================================================================================
*/

include { CONCAT_RNA } from '../modules/concat_rna_reads.nf'
include { CONCAT_DNA } from '../modules/concat_dna_reads.nf'

/*
========================================================================================
    Named workflow
========================================================================================
*/

workflow CONCATENATE {

if (params.process_rna){
    CONCAT_RNA(rna_reads)
}
if (params.process_dna){
    CONCAT_DNA(dna_reads)
}

}
