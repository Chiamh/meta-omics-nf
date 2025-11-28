//Concatenating R1 and R2 files that are spread across different lanes, using a common library ID
// This script assumes that the R1 and R2 files are in the same subdirectory, named after a LIBID like MHS584. e.g. $INPUTDIR/MHS584/*

process CONCAT_DNA {
	label "process_small"
	tag "${sample_id}"
	publishDir "${params.outdir}/DNA_merged"
	
	input:
	tuple val(sample_id), path(reads)
	
	output:
	tuple val(sample_id), path("${sample_id}_merged_1.fastq.gz"), path("${sample_id}_merged_2.fastq.gz")
	
	script:
	def r1_files = reads.findAll { it.name.matches('.*(_R1|_1)\\.(fq|fastq)(\\.gz)?$') }
	def r2_files = reads.findAll { it.name.matches('.*(_R2|_2)\\.(fq|fastq)(\\.gz)?$') }
	
	"""
	echo "Concatenating reads from ${sample_id}"
	
	cat ${r1_files.join(' ')} > "${sample_id}_merged_1.fastq.gz"
	
	cat ${r2_files.join(' ')} > "${sample_id}_merged_2.fastq.gz"
	echo "Finished processing ${sample_id}"
	"""


}
