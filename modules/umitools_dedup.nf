// Add UMI to header from RNAseq data using umitools

process UMITOOLS_DEDUP {
	label "process_medium"
	tag "${sample_id}"
	publishDir "${params.outdir}/decont/RNA", mode: 'copy'
	
	input:
	tuple val(sample_id), path(bam)
	val tag
	
	output:
	tuple val(sample_id), path("${sample_id}_${tag}.umidedup.bam"), emit: dedup_bam
	tuple val(sample_id), path("${sample_id}*_{1,2}.fastq.gz"), emit: dedup_fastq
	tuple val(sample_id), path("*.log"), emit: logs
	
	when:
	params.dedupe && params.process_rna
	
	script:
	"""
	umi_tools dedup -I ${bam}_bt2_pangenome_aligned.sorted.bam \\
	--output-stats=${sample_id}_${tag} -S ${sample_id}_${tag}.umidedup.bam
	
	samtools fastq -1 ${sample_id}_dedup_1.fastq.gz -2 ${sample_id}_dedup_2.fastq.gz ${sample_id}_${tag}.umidedup.bam
	"""
}

	
