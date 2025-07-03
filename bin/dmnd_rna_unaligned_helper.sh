#!/bin/bash
set -e
set -u

sample_id=${1}

grep -e ">" ${sample_id}_uniref90_unaligned.fa | awk 'sub(/^>/, "")' | awk '{ print $1 }' | sort | uniq > read_names
seqtk subseq ${sample_id}_mRNA_1.fastq.gz read_names > ${sample_id}_mRNA_1.unaligned.fastq && gzip ${sample_id}_mRNA_1.unaligned.fastq
seqtk subseq ${sample_id}_mRNA_2.fastq.gz read_names > ${sample_id}_mRNA_2.unaligned.fastq && gzip ${sample_id}_mRNA_2.unaligned.fastq
sed 's/ /_/g' ${sample_id}_uniref90_unaligned.fa > ${sample_id}_for_alignment.fa
mkdir ${sample_id}_unaligned && bowtie2-build ${sample_id}_for_alignment.fa ${sample_id}_unaligned/${sample_id}_unaligned
zcat ${sample_id}_mRNA_1.fastq.gz ${sample_id}_mRNA_2.fastq.gz | bowtie2 -q -x ${sample_id}_unaligned/${sample_id}_unaligned - | samtools view -bS | samtools sort > ${sample_id}.unaligned.bam && samtools index  ${sample_id}.unaligned.bam
