// Merge deduped PANGENOME- and DMND-aligned fastq for Kraken2 classification

process PAN_DMND_KRAKEN_FASTQ {
	label "process_medium"
	tag "${sample_id}"
	publishDir "${params.outdir}/kraken2_out/RNA", mode: 'copy'
	
	input:
	tuple val(sample_id), path(pan_umi_aligned_bam)
        tuple val(sample_id), path(dmnd_umi_aligned_bam)
        tuple val(sample_id), path(unaligned_umi_aligned_bam)
        tuple val(sample_id), path(dmnd_umi_aligned_out)
	tuple val(sample_id), path(reads)
	
	output:
	tuple val(sample_id), path("${sample_id}_pan_dmnd_merged_{1,2}.fastq.gz"), emit: merged
        tuple val(sample_id), path("${sample_id}_uniref90_aligned.out"), emit: dmnd        
	
	when:
	!params.profilers_off && params.process_rna

	script:
	"""
        mv MHS434_uniref90_aligned.out tmp
        samtools view ${dmnd_umi_aligned_bam} | cut -f1 | sort | uniq > read_names
        grep -f read_names tmp > MHS434_uniref90_aligned.out
        
        samtools view ${pan_umi_aligned_bam} | cut -f1 | sort | uniq >> read_names
        samtools view ${unaligned_umi_aligned_bam} | cut -f1 | sort | uniq >> read_names
        seqtk subseq ${reads[0]} read_names > ${sample_id}_pan_dmnd_merged_1.fastq && gzip ${sample_id}_pan_dmnd_merged_1.fastq
        seqtk subseq ${reads[1]} read_names > ${sample_id}_pan_dmnd_merged_2.fastq && gzip ${sample_id}_pan_dmnd_merged_2.fastq
	"""
}
