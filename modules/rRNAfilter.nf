// rRNA removal (human, fungal and bacteria) from RNAseq data using bbduk.sh
// If params.dedupe = false, then rRNA removal is the "final" step of decontamination
//The publishDir should not be included within a condition statement. Use instead a conditional expression like so:
// a ? b: c means if (a) b else c (ternary if special operator)

process RIBOFILTER {
	label "process_medium"
	label "error_retry"
	tag "${sample_id}"
	

	publishDir "${ params.dedupe ? params.outdir/decont/RNA/rRNAfilter_temp_fastq : params.outdir/decont/RNA }", mode: 'copy', pattern: '*.fastq.gz'
	publishDir "${params.outdir}/decont/RNA", mode: 'copy', pattern: '*.log'
	
	input:
	path ribokmers
	tuple val(sample_id), path(reads)
	
	output:
	tuple val(sample_id), path("${sample_id}*_{1,2}.fastq.gz"), emit: reads
	tuple val(sample_id), path("*.rRNAfilter.log") , emit: logs
	
	when:
	!params.decont_off && params.process_rna
	
	script:
	if (params.dedupe) {
	
	"""
	bbduk.sh in=${reads[0]} in2=${reads[1]} \\
	out=${sample_id}_mRNA_1.fastq.gz out2=${sample_id}_mRNA_2.fastq.gz \\
	k=31 \\
	ref=${ribokmers} stats=${sample_id}_rRNAfilter.log
	"""
	} else if (params.dedupe == false){
	"""
	bbduk.sh in=${reads[0]} in2=${reads[1]} \\
	out=${sample_id}_decont_1.fastq.gz out2=${sample_id}_decont_2.fastq.gz \\
	k=31 \\
	ref=${ribokmers} stats=${sample_id}_rRNAfilter.log
	"""
	}

}

	
	
