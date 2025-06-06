// Read de-duplication of RNAseq data using clumpify.sh

process DEDUP_CLUMPIFY {
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
	clumpify.sh in=${reads[0]} in2=${reads[1]} \\
	out=${sample_id}_decont_1.fastq.gz out2=${sample_id}_decont_2.fastq.gz \\
	dedupe=t \\
	optical=f 2>${sample_id}_dedup.log
	"""

}

	
