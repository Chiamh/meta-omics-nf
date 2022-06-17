//Concatenating R1 and R2 files that are spread across different lanes, using a common library ID
// This script assumes that the R1 and R2 files are in the same subdirectory, named after a LIBID like MHS584. e.g. $INPUTDIR/MHS584/*

process CONCAT_DNA {
	label "process_small"
	tag "${sample_id}"
	publishDir "${params.outdir}/DNA_merged"
	
	input:
	tuple val(sample_id), path(read1_files)
	tuple val(sample_id), path(read2_files)
	
	output:
	path("${sample_id}_merged_{1,2}.fastq.gz")
	
	script:
	"""
	echo concatenating reads from "${sample_id}"
	
	cat "${sample_id}"*_1.fq.gz > "${sample_id}_merged_1.fastq.gz"
	
	cat "${sample_id}"*_2.fq.gz > "${sample_id}_merged_2.fastq.gz"
	
	echo finished processing "${sample_id}"
	
	"""

}
