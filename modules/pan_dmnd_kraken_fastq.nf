// Merge deduped PANGENOME- and DMND-aligned fastq for Kraken2 classification

process PAN_DMND_KRAKEN_FASTQ {
	label "process_medium"
	tag "${sample_id}"
	publishDir "${params.outdir}/UMI_dedup_out/RNA", mode: 'copy'
	
	input:
	tuple val(sample_id), path(pan_umi_aligned_bam), path(dmnd_umi_aligned_bam), path(unaligned_dedup_readnames), path(dmnd_umi_aligned_out), path(reads)
	
	output:
	tuple val(sample_id), path("${sample_id}_pan_dmnd_merged_{1,2}.fastq.gz"), emit: merged
	tuple val(sample_id), path("${sample_id}_uniref90_aligned.out"), emit: dmnd
	tuple val(sample_id), path("${sample_id}_all_dedup_read_pair_stats.txt"), emit: logs
	
	when:
	!params.profilers_off && params.process_rna

	script:
	"""
        mv ${sample_id}_uniref90_aligned.out ${sample_id}_tmp
        samtools view ${dmnd_umi_aligned_bam} | cut -f1 | sort | uniq > ${sample_id}_dmnd_umi_read_names
        grep -Ff ${sample_id}_dmnd_umi_read_names ${sample_id}_tmp > ${sample_id}_uniref90_aligned.out
		
		wc -l ${sample_id}_dmnd_umi_read_names > ${sample_id}_all_dedup_read_pair_stats.txt
		
        samtools view ${pan_umi_aligned_bam} | cut -f1 | sort | uniq > ${sample_id}_pan_umi_read_names
		
		wc -l ${sample_id}_pan_umi_read_names >> ${sample_id}_all_dedup_read_pair_stats.txt
		
		cat ${sample_id}_pan_umi_read_names ${sample_id}_dmnd_umi_read_names | sort | uniq > ${sample_id}_all_aligned_uniq_read_names
		
		wc -l ${sample_id}_all_aligned_uniq_read_names >> ${sample_id}_all_dedup_read_pair_stats.txt
		
        cat ${unaligned_dedup_readnames} ${sample_id}_all_aligned_uniq_read_names | sort | uniq > ${sample_id}_read_pairs_to_extract
		
        seqtk subseq ${reads[0]} ${sample_id}_read_pairs_to_extract > ${sample_id}_pan_dmnd_merged_1.fastq && gzip ${sample_id}_pan_dmnd_merged_1.fastq
        seqtk subseq ${reads[1]} ${sample_id}_read_pairs_to_extract > ${sample_id}_pan_dmnd_merged_2.fastq && gzip ${sample_id}_pan_dmnd_merged_2.fastq
		
		rm ${sample_id}_tmp
	"""
}
