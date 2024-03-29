// Bracken classification for metagenomic data. Should not be run on metatranscriptomes

/*
params.kraken2db is the path to a built Kraken2 database which also contains the Bracken database files
*/

process BRACKEN {
	label "process_small"
	tag "${sample_id}"
	publishDir "${params.outdir}/kraken2_out/DNA", mode: 'copy'
	
	input:
	path kraken2db
	val(readlength)
	tuple val(sample_id), path("${sample_id}_kraken2.tax")
	
	output:
	tuple val(sample_id), path("${sample_id}_bracken.tax"), emit: results
	tuple val(sample_id), path("${sample_id}_bracken.out"), emit: out

	when:
	!params.profilers_off && params.process_dna
	
	script:
	"""	
	bracken -d "${kraken2db}" -i "${sample_id}_kraken2.tax" \\
	-o "${sample_id}_bracken.out" -w "${sample_id}_bracken.tax" -r "${readlength}" -l S
		
	"""
}
