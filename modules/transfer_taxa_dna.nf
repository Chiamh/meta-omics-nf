// Transfer kraken2 taxonomic labels for functional annotation of MGX data
process TRF_TAXA_DNA {
	label "process_medium"
	tag "${sample_id}"
	publishDir "${params.outdir}/MGX_annotations", mode: 'copy'
	
	input:
	path pangenome_annots
	tuple val(sample_id), path(k2_out), path(panalign_bam), path(dmnd_aligned_results), path(unaligned_fa)
	
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
 
