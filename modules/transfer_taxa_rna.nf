// Transfer kraken2 taxonomic labels for functional annotation of MTX data
process TRF_TAXA_RNA {
	label "process_medium"
	tag "${sample_id}"
	publishDir "${params.outdir}/MTX_annotations", mode: 'copy'
	
	input:
	path pangenome_annots
	tuple val(sample_id), path(k2_out)
	tuple val(sample_id), path(panalign_bam)
	tuple val(sample_id), path(dmnd_aligned_results)
	tuple val(sample_id), path(unaligned_fa)
	
	output:
	tuple val(sample_id), path("${sample_id}_all_aligned_taxonomy.tsv"), emit: aligned_taxa
	tuple val(sample_id), path("${sample_id}_unaligned_taxonomy.tsv"), emit: unaligned_taxa
	
	when:
	!params.annotate_off
	
	script:

	"""
	
	transfer_taxa_helper.sh "${sample_id}" "${pangenome_annots}"

	"""
}
 
