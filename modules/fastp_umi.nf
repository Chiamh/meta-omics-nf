// Add UMI to header from RNAseq data using fastp

process FASTP_UMI {
	label "process_medium"
	tag "${sample_id}"
	publishDir "${params.outdir}/decont/RNA", mode: 'copy'
	
	input:
	tuple val(sample_id), path(reads_file)
	
	output:
	tuple val(sample_id), path("${sample_id}_fastp_{1,2}.fastq.gz"), emit: reads
	tuple path("${sample_id}.html"), path("${sample_id}.json") , emit: logs
	
	when:
	params.process_rna
	
	script:
	"""
	
	fastp -i ${reads_file[0]} -I ${reads_file[1]} \\
	--out1 ${sample_id}_extract_1.fastq.gz --out2 ${sample_id}_extract_2.fastq.gz \\
	-U --umi_loc=read1 --umi_len=11 --umi_prefix UMI -j ${sample_id}.json -h ${sample_id}.html
	
	zcat ${sample_id}_extract_1.fastq.gz | sed 's/:UMI_/_/g' | gzip > ${sample_id}_fastp_1.fastq.gz
	
	rm ${sample_id}_extract_1.fastq.gz
	
	zcat ${sample_id}_extract_2.fastq.gz | sed 's/:UMI_/_/g' | gzip > ${sample_id}_fastp_2.fastq.gz
	
	rm ${sample_id}_extract_2.fastq.gz
	
	"""
}

	
