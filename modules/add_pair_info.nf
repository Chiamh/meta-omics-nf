// Add UMI to header from RNAseq data using umitools

process ADD_PAIR_INFO {
	label "process_medium"
	tag "${sample_id}"
	publishDir "${params.outdir}/decont/RNA", mode: 'copy'
	
	input:
	tuple val(sample_id), path(reads)
	
	output:
	tuple val(sample_id), path("${sample_id}*_{1,2}.fastq.gz"), emit: reads
	
	when:
	params.dedupe && params.process_rna
	
	script:
        """
        add_pair_info.sh ${reads[0]} ${reads[1]}
        """
}

	
