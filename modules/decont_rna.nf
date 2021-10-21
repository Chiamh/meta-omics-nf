// decontamination/removal of human reads using STAR aligner, followed by rRNA removal and deduplication to get microbial reads

process DECONT_RNA {
	label "process_high"
	tag "${sample_id}"
	publishDir "${params.outdir}/decont/RNA", mode: 'copy'
	
	input:
	path star_index
	path ribokmers
	tuple val(sample_id), path(reads_file)
	
	output:
	tuple val(sample_id), path("*.fastq.gz"), emit: reads
	tuple path("${sample_id}.html"), path("${sample_id}.json"), path("${sample_id}_rRNAfilter.log"), path("${sample_id}_dedup.log") , emit: logs
	
	when:
	!params.decont_off && params.process_rna
	
	script:
	
	if (params.dedupe) {
		"""
		fastp -i ${reads_file[0]} -I ${reads_file[1]} \\
		--out1 ${sample_id}_fastp_1.fastq.gz --out2 ${sample_id}_fastp_2.fastq.gz \\
		-j ${sample_id}.json -h ${sample_id}.html
		
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
		
		rm ${sample_id}_fastp_1.fastq.gz
		rm ${sample_id}_fastp_2.fastq.gz
	
		bbduk.sh in=${sample_id}_unmapped_1.fastq.gz in2=${sample_id}_unmapped_2.fastq.gz \\
		out=${sample_id}_mRNA_1.fastq.gz out2=${sample_id}_mRNA_2.fastq.gz \\
		k=31 \\
		ref=${ribokmers} \\
		stats=${sample_id}_rRNAfilter.log
		
		rm ${sample_id}_unmapped_1.fastq.gz
		rm ${sample_id}_unmapped_2.fastq.gz
	
		clumpify.sh in=${sample_id}_mRNA_1.fastq.gz in2=${sample_id}_mRNA_2.fastq.gz \\
		out=${sample_id}_decont_1.fastq.gz out2=${sample_id}_decont_2.fastq.gz \\
		dedupe=t \\
		optical=f 2>${sample_id}_dedup.log
	
		rm ${sample_id}_mRNA_1.fastq.gz
		rm ${sample_id}_mRNA_2.fastq.gz
		"""
		} else if (!params.dedupe){
		"""
		fastp -i ${reads_file[0]} -I ${reads_file[1]} \\
		--out1 ${sample_id}_fastp_1.fastq.gz --out2 ${sample_id}_fastp_2.fastq.gz \\
		-j ${sample_id}.json -h ${sample_id}.html
		
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
		
		rm ${sample_id}_fastp_1.fastq.gz
		rm ${sample_id}_fastp_2.fastq.gz
	
		bbduk.sh in=${sample_id}_unmapped_1.fastq.gz in2=${sample_id}_unmapped_2.fastq.gz \\
		out=${sample_id}_decont_1.fastq.gz out2=${sample_id}_decont_2.fastq.gz \\
		k=31 \\
		ref=${ribokmers} \\
		stats=${sample_id}_rRNAfilter.log
		
		rm ${sample_id}_unmapped_1.fastq.gz
		rm ${sample_id}_unmapped_2.fastq.gz
	
		"""
		}
}

