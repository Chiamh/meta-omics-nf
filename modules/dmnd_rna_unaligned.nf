// Translated alignment search for metatranscriptomes

/*
Requires a pre-built diamond2 database. Can use uniref90 or similar. 
*/

process DMND_RNA_UNALIGNED {
	label "process_highmem"
	tag "${sample_id}"
	publishDir "${params.outdir}/MTX_dmnd_out", mode: 'copy'
	
	input:
	tuple val(sample_id), path(reads)
        tuple val(sample_id), path(unaligned_fastq)
	
	output:
	tuple val(sample_id), path("${sample_id}.unaligned.bam"), path("${sample_id}.unaligned.bam.bai"), emit: bam

	when:
	!params.panalign_off && !params.diamond_off
	
	script:
	"""	
        dmnd_rna_unaligned_helper.sh "${sample_id}"
	"""
}
