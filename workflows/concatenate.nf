// WORKFLOW file for concatenating R1 and R2 fastq.gz files that are spread across different lanes, using a common library ID
// Assumes that the R1 and R2 files are in the same subdirectory, named after a LIBID like MHS584. e.g. $INPUTDIR/MHS584/*
// params.rna_reads and params.dna_reads are paths to folders containing these subdirectories
// Yes you could create a glob pattern parameter to make this more generalizable, but this meets my needs for now.
/*
========================================================================================
    Define channels for read pairs
========================================================================================
*/

if (params.process_rna){
    libid_rna = Channel.fromPath( [params.rna_reads + '/**MHS[0-9]*'], checkIfExists:true ).map{ file ->
	libid = file.getParent().getName(); return libid //the parent dir is likely the library ID
	}
	.unique()
}

if (params.process_dna){
    libid_dna = Channel.fromPath( [params.dna_reads + '/**MHS[0-9]*'], checkIfExists:true ).map{ file ->
	libid = file.getParent().getName(); return libid //the parent dir is likely the library ID
	}
	.unique()
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

workflow CONCAT {

if (params.process_rna){
    CONCAT_RNA(params.rna_reads, libid_rna)
}
if (params.process_dna){
    CONCAT_DNA(params.dna_reads, libid_dna)
}

}