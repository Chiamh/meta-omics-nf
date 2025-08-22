// Annotation of pangenome aligned clusters (e.g. IHSMGC) which have sufficient MGX read coverage, with Uniref90 labels, Eggnog COG labels and GO terms 
// No need to run eggnog mapper here because all the pangenes have been annotated by eggnog already (pangenome_annots)
process ANNOT_PAN_DNA {
	label "process_high"
	tag "${sample_id}"
	publishDir "${params.outdir}/MGX_annotations", mode: 'copy'
	
	input:
	path pangenome_annots
	tuple val(sample_id), path(panalign_results)
	
	output:
	tuple val(sample_id), path("${sample_id}_panalign_annot.tsv"), emit: results
	
	when:
	!params.annotate_off
	
	script:

	"""
	
	annot_pan_helper.sh "${sample_id}" "${pangenome_annots}"
	
	"""
}
 
