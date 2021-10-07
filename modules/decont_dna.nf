// decontamination or removal of human reads from metagenomes

params.hg_fasta = './genomes/hg38_with_IVT/hg38_with_IVT.fa'

process DECONT_DNA {
	label "process_medium"
	label "error_retry"
	tag "${sample_id}"
	publishDir "${params.outdir}/decont/DNA", mode: 'copy'
	
	input:
	path hg_fasta
	tuple val(sample_id), path(reads_file)
	
	output:
	tuple val(sample_id), path("*.fastq.gz"), emit: reads
	tuple path("${sample_id}.html"), path("${sample_id}.json") , emit: logs
	
	when:
	!params.decont_off && params.process_dna
	
	script:
	"""
	fastp -i ${reads_file[0]} -I ${reads_file[1]} --stdout -j ${sample_id}.json -h ${sample_id}.html | \
	bwa mem -p -t $task.cpus ${hg_fasta} - | \
	samtools fastq -f12 -F256 -1 ${sample_id}_decont_1.fastq.gz -2 ${sample_id}_decont_2.fastq.gz -
		
	"""
}

