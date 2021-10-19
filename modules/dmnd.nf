// Translated alignment search for metatranscriptomes

/*
Requires a pre-built diamond2 database. Can use uniref90 or similar. 
*/

process DMND {
	label "process_high"
	label "error_retry"
	tag "${sample_id}"
	publishDir "${params.outdir}/dmnd_out", mode: 'copy'
	
	input:
	path dmnddb
	tuple val(sample_id), path(read_file)
	
	output:
	tuple val(sample_id), path("${sample_id}_uniref90_aligned.out"), emit: aligned
	tuple val(sample_id), path("${sample_id}_uniref90_unaligned.fa"), emit: unaligned
	tuple val(sample_id), path("${sample_id}_dmnd.log"), emit: logs

	when:
	!params.panalign_off && !params.diamond_off
	
	script:
	"""	
	diamond blastx --query "${read_file}" \
	--id 80 --query-cover 90 --threads $task.cpus --max-target-seqs 1 \
	-b6 \
	--outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen evalue bitscore salltitles \
	--db "${dmnddb}" \
	--out "${sample_id}_uniref90_aligned.out" \
	--un "${sample_id}_uniref90_unaligned.fa" \
	--unfmt fasta 2>"${sample_id}_dmnd.log"
	"""
}