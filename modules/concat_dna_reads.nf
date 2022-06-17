//Concatenating R1 and R2 files that are spread across different lanes, using a common library ID
// This script assumes that the R1 and R2 files are in the same subdirectory, named after a LIBID like MHS584. e.g. $INPUTDIR/MHS584/*

process CONCAT_DNA {
	label "process_small"
	tag "${sample_id}"
	publishDir "${params.outdir}/DNA_merged", mode: 'move'
	
	input:
	path readsdir
	val sample_id
	
	output:
	path("${sample_id}_merged_{1,2}.fastq.gz")
	
	script:
	"""
	echo concatenating reads from "${sample_id}"
	
	cat ${readsdir}/${sample_id}/*_1.fq.gz > "${sample_id}_merged_1.fastq.gz"
	
	cat ${readsdir}/${sample_id}/*_2.fq.gz > "${sample_id}_merged_2.fastq.gz"
	
	echo finished processing "${sample_id}"
	
	"""

}