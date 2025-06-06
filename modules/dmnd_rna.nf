// Translated alignment search for metatranscriptomes

/*
Requires a pre-built diamond2 database. Can use uniref90 or similar. 
*/

process DMND_RNA {
	label "process_highmem"
	tag "${sample_id}"
	publishDir "${params.outdir}/MTX_dmnd_out", mode: 'copy'
	
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
	diamond blastx --query "${read_file}" \\
	--id 80 --query-cover 90 --threads $task.cpus --max-target-seqs 1 \\
	-b6 \\
	--outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen evalue bitscore salltitles \\
	--db "${dmnddb}" \\
	--out "${sample_id}_uniref90_aligned.out" \\
	--un "${sample_id}_uniref90_unaligned.fa" \\
	--unfmt fasta 2>"${sample_id}_dmnd.log"

	grep "UniRef90" ${sample_id}_uniref90_aligned.out | awk '{print $3}' > ${sample_id}_uniref90ids.txt
	head -n 1 ${sample_id}_uniref90_aligned.out > ${sample_id}.sam
	grep -Ff ${sample_id}_uniref90ids.txt ${db_dir}/uniref90_09Jun2021/uniref90.fai | awk '{print "@SQ\tSN:" $1 "\tLN:" $2}' >> ${sample_id}.sam
	tail -n +2 ${sample_id}_uniref90_aligned.out >> ${sample_id}.sam
	sed -i 's/PN:DIAMOND/ID:DIAMOND/' ${sample_id}.sam
	sed -i '/@mm/d' ${sample_id}.sam

	samtools sort ${sample_id}.sam > ${sample_id}.sorted.bam
	samtools index ${sample_id}.sorted.bam
	"""
}
