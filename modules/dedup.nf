// Read de-duplication of RNAseq data using clumpify.sh

process DEDUP {
	label "process_medium"
	label "error_retry"
	tag "${sample_id}"
	publishDir "${params.outdir}/decont/RNA", mode: 'copy'
	
	input:
	tuple val(sample_id), path(reads_file)
	
	output:
	tuple val(sample_id), path("*.fastq.gz"), emit: reads
	tuple val(sample_id), path("*.log"), emit: logs
	
	when:
	!params.decont_off && params.dedupe && params.process_rna
	
	script:
	"""
	clumpify.sh in=${reads_file[0]} in2=${reads_file[1]} \
	out=${sample_id}_decont_1.fastq.gz out2=${sample_id}_decont_2.fastq.gz \
	dedupe=t \
	optical=f 2>${sample_id}_dedup.log
	"""

}

	