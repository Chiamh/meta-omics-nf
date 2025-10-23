// De-duplication using UMItools API: UMIclusterer

process UMITOOLS_CLUST {
	label "process_medium"
	tag "${sample_id}"
	publishDir "${params.outdir}/decont/RNA", mode: 'copy'
	
	input:
	tuple val(sample_id), path(clusters)
	val tag
	
	output:
	tuple val(sample_id), path("${sample_id}_umiclusterer_output.tsv"), emit: umiclustout
	tuple val(sample_id), path("${sample_id}_read_names_from_unaligned_dedup"), emit: readnames
	
	when:
	params.dedupe && params.process_rna
	
	script:
	"""
	bash format_clstr_file_after_vsearch.sh ${sample_id}_uniref90_unaligned_vsearch_0.95.out > ${sample_id}_barcode_by_cluster.txt
	
	bash count_clstrs.sh ${sample_id}_barcode_by_cluster.txt > ${sample_id}_input_to_umiclusterer.tsv
	
	python run_umi_clusterer.py ${sample_id}_input_to_umiclusterer.tsv ${sample_id}_umiclusterer_output.tsv directional
	
	bash extract_deduped_read_IDs_with_umi.sh ${sample_id}_umiclusterer_output.tsv <(awk -F "\t" -v OFS="\t" '{print \$9, \$2}' ${sample_id}_uniref90_unaligned_vsearch_0.95.out ) ${sample_id}
	
	"""
}

