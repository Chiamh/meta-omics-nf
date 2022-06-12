// Annotation of candidate uniref90 clusters which have sufficient MTX read coverage, with Eggnog COG labels and GO terms 

process ANNOT_DMND_RNA {
	label "process_high"
	tag "${sample_id}"
	publishDir "${params.outdir}/MTX_annotations", mode: 'copy'
	
	input:
	path uniref90_fasta
	path eggnog_OG_annots
	path eggnog_db
	path uniref90_GO
	tuple val(sample_id), path(dmnd_results)
	
	output:
	tuple val(sample_id), path("${sample_id}.emapper.annotations"), path("${sample_id}_transl-search_annot.tsv"), emit: results
	
	when:
	!params.annotate_off
	
	script:

	"""
	annot_dmnd_helper_1.sh "${sample_id}" "${dmnd_results}" "${uniref90_fasta}"
	
	emapper.py --cpu ${task.cpus} -i "${sample_id}"_chosen.fa --data_dir "${eggnog_db}" -m diamond --go_evidence all --output "${sample_id}" --output_dir ./ --dbmem

	annot_dmnd_helper_2.sh "${sample_id}"
	
	extract_OG.sh "${sample_id}"_eggnog_index "${sample_id}"_eggnog_OGs "${sample_id}"
	
	annot_dmnd_helper_3.sh "${sample_id}" "${eggnog_OG_annots}" "${dmnd_results}" "${uniref90_GO}"
	
	"""

}
	
	




