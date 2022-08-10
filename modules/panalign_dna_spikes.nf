// Pangenome/gene catalog nucleotide alignment for metagenomes (MGX) with microbial spike in removal

/*
Requires a pre-built pangenome with a bowtie2 index
This process will concatenate the metagenomic R1 and R2 files before bowtie2 mapping to pangenome.
*/

process PANALIGN_DNA_SPIKES {
	label "process_highmem"
	tag "${sample_id}"
	publishDir "${params.outdir}/MGX_panalign_out", mode: 'copy'
	
	input:
	path pangenome_path
	path spike_in_path
	tuple val(sample_id), path(reads_file)
	tuple val(sample_id), path(k2_out)
	
	output:
	tuple val(sample_id), path("${sample_id}_bt2_pangenome_aligned.bam"), emit: aligned
	tuple val(sample_id), path("${sample_id}_bt2_pangenome_aligned_filtered_cov.tsv"), emit: coverage
	tuple val(sample_id), path("${sample_id}_bt2_pangenome_unaligned.fastq.gz"), emit: unaligned
	tuple val(sample_id), path("${sample_id}_bt2.log"), emit: logs
	
	when:
	!params.panalign_off && params.process_dna && params.rm_spikes
	
	script:
	"""
	zcat ${reads_file[0]} ${reads_file[1]} | \\
	(bowtie2 -q -x ${pangenome_path}/${params.pangenome} -U - --un-gz "${sample_id}_bt2_pangenome_unaligned_unfiltered.fastq.gz" \\
	-p $task.cpus --very-sensitive) 2>"${sample_id}_bt2.log" | \\
	samtools view -bS - > "${sample_id}_bt2_pangenome_aligned_unfiltered.bam"

	panalign_spikes_helper.sh "${sample_id}" "${spike_in_path}"

		
	"""
}
