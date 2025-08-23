// Add UMI to header from RNAseq data using umitools

process UMITOOLS_DEDUP {
	label "process_high"
	tag "${sample_id}"
	publishDir "${params.outdir}/decont/RNA", mode: 'copy'
	
	input:
	tuple val(sample_id), path(bam), path(bai)
	val tag
	
	output:
	tuple val(sample_id), path("${sample_id}_${tag}.umidedup.bam"), emit: dedup_bam
	tuple val(sample_id), path("*.log"), emit: logs
	
	when:
	params.dedupe && params.process_rna
	
	script:
	"""
	umi_tools dedup -I ${bam} --paired \\
	--output-stats=${sample_id}_${tag} -S ${sample_id}_${tag}.umidedup.bam > ${sample_id}_${tag}_umi_dedup.log
	"""
}

	
