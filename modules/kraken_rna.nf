// Kraken2 classification

/*
It is encouraged for users to run kraken2 with 10-20 threads. 
*/

process KRAKEN2_RNA {
	label "process_high"
	tag "${sample_id}"
	publishDir "${params.outdir}/kraken2_out/RNA", mode: 'copy'
	
	input:
	path kraken2db
	tuple val(sample_id), path(reads_file)
	
	output:
	tuple val(sample_id), path("${sample_id}_kraken2.tax"), emit: k2tax
	tuple val(sample_id), path("${sample_id}_kraken2.out"), emit: k2out
	tuple val(sample_id), path("${sample_id}_k2.s.tsv"), emit: speciesreport
	tuple val(sample_id), path("${sample_id}_kraken2_minimizer.tax"), emit: k2mintax
	
	when:
	!params.profilers_off && params.process_rna

	script:
	"""
	kraken2 --use-names --threads ${task.cpus} --db "${kraken2db}" \\
	--report "${sample_id}_kraken2_minimizer.tax" --report-minimizer-data --output "${sample_id}_kraken2.out" \\
	--gzip-compressed --paired ${reads_file[0]} ${reads_file[1]}
	
	k2_helper.sh "${sample_id}"
	
	"""
}
