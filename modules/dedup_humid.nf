// Read de-duplication of RNAseq data using humid

process DEDUP_HUMID {
    label "process_medium"
    tag "${sample_id}"
    publishDir "${params.outdir}/decont/RNA", mode: 'copy'
	
    input:
    tuple val(sample_id), path(reads)
	
    output:
    tuple val(sample_id), path("${sample_id}*_{1,2}_depup.fastq.gz"), emit: reads //there's also *_annotated.fastq.gz
    tuple val(sample_id), path("*.log"), emit: logs
	
    when:
    params.dedupe && params.process_rna
	
    script:
    """
    humid -d ${sample_id} -l ${sample_id}.log ${reads[0]} ${reads[1]}
    mv ${sample_id}/*.fastq* .
    """
}

	
