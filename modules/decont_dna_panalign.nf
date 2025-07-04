// decontamination or removal of human reads from metagenomes

process DECONT_DNA_PANALIGN {
	label "process_high"
	tag "${sample_id}"
	publishDir "${params.outdir}/decont/DNA", mode: 'copy'
	
	input:
	path human_pangenome_path
	tuple val(sample_id), path(reads_file)

	output:
	tuple val(sample_id), path("${sample_id}*_{1,2}.fastq.gz"), emit: reads
	tuple val(sample_id), path("${sample_id}_hprc-v1.1-mc-chm13_QC.txt"), emit: stats
	
	when:
	!params.decont_off && params.process_dna
	
	script:
	"""
        decont_dna_panalign_helper.sh "${sample_id}" "${reads_file[0]}" "${reads_file[1]}" "${human_pangenome_path}" "${task.cpus}"
	"""
}
