// decontamination or removal of human reads from metagenomes

process DECONT_DNA {
	label "process_high"
	label "error_retry"
	tag "${sample_id}"
	publishDir "${params.outdir}/decont/DNA", mode: 'copy'
	
	input:
	path bwaidx_path
	tuple val(sample_id), path(reads_file)
	
	output:
	tuple val(sample_id), path("*.fastq.gz"), emit: reads
	tuple path("${sample_id}.html"), path("${sample_id}.json") , emit: logs
	
	when:
	!params.decont_off && params.process_dna
	
	script:
	"""
	fastp -i ${reads_file[0]} -I ${reads_file[1]} --stdout -j ${sample_id}.json -h ${sample_id}.html | \\
	bwa mem -p -t $task.cpus ${bwaidx_path}/${params.bwaidx} - | \\
	samtools fastq -f12 -F256 -1 ${sample_id}_decont_1.fastq.gz -2 ${sample_id}_decont_2.fastq.gz -
		
	"""
}

