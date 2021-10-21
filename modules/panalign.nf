// Pangenome/gene catalog nucleotide alignment for metatranscriptomes

/*
Requires a pre-built pangenome with a bowtie2 index
This process will concatenate the metatranscriptomic R1 and R2 files before bowtie2 mapping to pangenome.
*/

process PANALIGN {
	label "process_medium"
	label "error_retry"
	tag "${sample_id}"
	publishDir "${params.outdir}/panalign_out", mode: 'copy'
	
	input:
	path pangenome_path
	tuple val(sample_id), path(read_file)
	
	output:
	tuple val(sample_id), path("${sample_id}.bt2_pangenome_aligned.bam"), emit: aligned
	tuple val(sample_id), path("${sample_id}.bt2_pangenome_unaligned.fastq.gz"), emit: unaligned
	tuple val(sample_id), path("${sample_id}_bt2.log"), emit: logs
	
	when:
	!params.panalign_off && params.process_rna
	
	script:
	"""	
	zcat read_file[0] read_file[1] | \\
	(bowtie2 -q -x ${pangenome_path}/${params.pangenome} -U - --un-gz "${sample_id}_bt2_pangenome_unaligned.fastq.gz" \\
	-p $task.cpus --very-sensitive) 2>"${sample_id}_bt2.log" | \\
	samtools view -bS - > "${sample_id}_bt2_pangenome_aligned.bam"; done
		
	"""
}
