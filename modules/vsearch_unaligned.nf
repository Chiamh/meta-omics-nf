// Vsearch clustering of reads

process VSEARCH {
	label "process_high"
	tag "${sample_id}"
	publishDir "${params.outdir}/MTX_dmnd_out", mode: 'copy'
	
	input:
	tuple val(sample_id), path(unaligned_fastq)
	
	output:
	tuple val(sample_id), path("${sample_id}_uniref90_unaligned_vsearch_0.95.out"), emit: clusters

	when:
	!params.panalign_off && !params.diamond_off
	
	script:
	"""	
	vsearch --cluster_fast ${sample_id}_uniref90_unaligned.fa --clusterout_id --id 0.95 --iddef 2 --uc ${sample_id}_uniref90_unaligned_vsearch_0.95.out
	"""
}
