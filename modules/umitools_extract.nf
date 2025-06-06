// Add UMI to header from RNAseq data using umitools

process UMITOOLS_EXTRACT {
	label "process_medium"
	tag "${sample_id}"
	publishDir "${params.outdir}/decont/RNA", mode: 'copy'
	
	input:
	tuple val(sample_id), path(reads)
	
	output:
	tuple val(sample_id), path("${sample_id}*_{1,2}.fastq.gz"), emit: reads
	tuple val(sample_id), path("*.log"), emit: logs
	
	when:
	params.dedupe && params.process_rna
	
	script:
	"""
	umi_tools extract --bc-pattern=NNNNNNNNNNN \\
	-I ${reads[0]} -S ${sample_id}_decont_1.fastq.gz \\
	--read2-in=${reads[1]} --read2-out=${sample_id}_decont_2.fastq.gz
	"""

}

	
