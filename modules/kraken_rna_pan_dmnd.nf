// Kraken2 classification

/*
It is encouraged for users to run kraken2 with 10-20 threads. 
*/

process KRAKEN2_RNA_PAN_DMND {
	label "process_high"
	tag "${sample_id}"
	publishDir "${params.outdir}/kraken2_out/RNA", mode: 'copy'
	
	input:
	path kraken2db
	tuple val(sample_id), path(pan_reads_file), path(dmnd_reads_file)
	
	output:
	tuple val(sample_id), path("${sample_id}_pan_kraken2.tax"), emit: pan_k2tax
	tuple val(sample_id), path("${sample_id}_dmnd_kraken2.tax"), emit: dmnd_k2tax
	tuple val(sample_id), path("${sample_id}_pan_kraken2.out"), emit: pan_k2out
	tuple val(sample_id), path("${sample_id}_dmnd_kraken2.out"), emit: dmnd_k2out
	tuple val(sample_id), path("${sample_id}_pan_k2.s.tsv"), emit: pan_speciesreport
	tuple val(sample_id), path("${sample_id}_dmnd_k2.s.tsv"), emit: dmnd_speciesreport
	tuple val(sample_id), path("${sample_id}_pan_kraken2_minimizer.tax"), emit: pan_k2mintax
	tuple val(sample_id), path("${sample_id}_dmnd_kraken2_minimizer.tax"), emit: dmnd_k2mintax
	tuple val(sample_id), path("${sample_id}_pan_dmnd_kraken2.out"), emit: k2out
	
	when:
	!params.profilers_off && params.process_rna

	script:
	"""
	kraken2 --use-names --threads ${task.cpus} --db "${kraken2db}" \\
	--report "${sample_id}_pan_kraken2_minimizer.tax" --report-minimizer-data --output "${sample_id}_pan_kraken2.out" \\
	--gzip-compressed --paired ${pan_reads_file[0]} ${pan_reads_file[1]}

	k2_helper.sh "${sample_id}_pan"

	kraken2 --use-names --threads ${task.cpus} --db "${kraken2db}" \\
	--report "${sample_id}_dmnd_kraken2_minimizer.tax" --report-minimizer-data --output "${sample_id}_dmnd_kraken2.out" \\
	--gzip-compressed --paired ${dmnd_reads_file[0]} ${dmnd_reads_file[1]}

	k2_helper.sh "${sample_id}_dmnd"

	cat ${sample_id}_pan_kraken2.out ${sample_id}_dmnd_kraken2.out > ${sample_id}_pan_dmnd_kraken2.out
	"""
}
