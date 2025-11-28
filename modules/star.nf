// decontamination or removal of human reads from RNAseq using STAR. 
// a ? b: c means if (a) b else c (ternary if special operator)  
process STAR {
	label "process_high"
	tag "${sample_id}"
	publishDir { !params.dedupe && !params.remove_rRNA ? "${params.outdir}/decont/RNA" : "${params.outdir}/decont/RNA/tmp_fastq" }, mode: 'copy', pattern: '*_unmapped_{1,2}.fastq.gz'
	
	
	input:
	path star_index
	tuple val(sample_id), path(reads_file)
	
	output:
	tuple val(sample_id), path("${sample_id}_unmapped_{1,2}.fastq.gz"), emit: microbereads
	
	when:
	!params.decont_off && params.process_rna
	
	script:
	"""
		
		STAR --runMode alignReads \\
			 --runThreadN $task.cpus \\
			 --outSAMtype None \\
			 --readFilesCommand zcat \\
			 --genomeDir ${star_index} \\
			 --outFileNamePrefix ${sample_id}. \\
			 --readFilesIn ${sample_id}_fastp_1.fastq.gz ${sample_id}_fastp_2.fastq.gz \\
			 --outReadsUnmapped Fastx
		
		if [ -f ${sample_id}.Unmapped.out.mate1 ]; then
        mv ${sample_id}.Unmapped.out.mate1 ${sample_id}_unmapped_1.fastq
        gzip ${sample_id}_unmapped_1.fastq
		fi
    	if [ -f ${sample_id}.Unmapped.out.mate2 ]; then
        mv ${sample_id}.Unmapped.out.mate2 ${sample_id}_unmapped_2.fastq
        gzip ${sample_id}_unmapped_2.fastq
    	fi
		
	"""
}

