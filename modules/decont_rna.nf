// decontamination/removal of human reads, followed by rRNA removal and deduplication

params.hg_fasta = './genomes/hg38_with_IVT/hg38_with_IVT.fa'
params.ribokmers = './databases/ribokmers.fa.gz'

process DECONT_RNA {
	label "process_medium"
	tag "${sample_id}"
	publishDir "${params.outdir}/decont/RNA", mode: 'copy'
	
	input:
	path hg_fasta
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
		fastp -i ${reads_file[0]} -I ${reads_file[1]} --stdout -j ${sample_id}.json -h ${sample_id}.html | \
		bwa mem -p -t $task.cpus ${hg_fasta} - | \
		samtools fastq -f12 -F256 -1 ${sample_id}_fastp_1.fastq.gz -2 ${sample_id}_fastp_2.fastq.gz -
	
		bbduk.sh in=${sample_id}_fastp_1.fastq.gz in2=${sample_id}_fastp_1.fastq.gz \
		out=${sample_id}_mRNA_1.fastq.gz out2=${sample_id}_mRNA_2.fastq.gz \
		k=31 \
		ref=${ribokmers} \
		stats=${sample_id}_rRNAfilter.log
	
		rm ${sample_id}_fastp_1.fastq.gz
		rm ${sample_id}_fastp_2.fastq.gz
	
		clumpify.sh in=${sample_id}_mRNA_1.fastq.gz in2=${sample_id}_mRNA_2.fastq.gz \
		out=${sample_id}_decont_1.fastq.gz out2=${sample_id}_decont_2.fastq.gz \
		dedupe=t \
		optical=f 2>${sample_id}_dedup.log
	
		rm ${sample_id}_mRNA_1.fastq.gz
		rm ${sample_id}_mRNA_2.fastq.gz
		"""
		} else if (!params.dedupe){
		"""
		fastp -i ${reads_file[0]} -I ${reads_file[1]} --stdout -j ${sample_id}.json -h ${sample_id}.html | \
		bwa mem -p -t $task.cpus ${hg_fasta} - | \
		samtools fastq -f12 -F256 -1 ${sample_id}_fastp_1.fastq.gz -2 ${sample_id}_fastp_2.fastq.gz -
	
		bbduk.sh in=${sample_id}_fastp_1.fastq.gz in2=${sample_id}_fastp_1.fastq.gz \
		out=${sample_id}_decont_1.fastq.gz out2=${sample_id}_decont_2.fastq.gz \
		k=31 \
		ref=${ribokmers} \
		stats=${sample_id}_rRNAfilter.log
	
		rm ${sample_id}_fastp_1.fastq.gz
		rm ${sample_id}_fastp_2.fastq.gz	
		"""
		}
}

