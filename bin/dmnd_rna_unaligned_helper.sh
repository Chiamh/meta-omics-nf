#!/bin/bash
set -e
set -u

sample_id=${1}

grep -e ">" ${sample_id}_uniref90_unaligned.fa | awk 'sub(/^>/, "")' | awk '{ print $1 }' | sort | uniq > read_names
seqtk subseq ${sample_id}_mRNA_1.fastq.gz read_names > ${sample_id}_mRNA_1.unaligned.fastq && gzip ${sample_id}_mRNA_1.unaligned.fastq
seqtk subseq ${sample_id}_mRNA_2.fastq.gz read_names > ${sample_id}_mRNA_2.unaligned.fastq && gzip ${sample_id}_mRNA_2.unaligned.fastq
samtools import -1 "${sample_id}"_mRNA_1.unaligned.fastq.gz -2 "${sample_id}"_mRNA_2.unaligned.fastq.gz -O sam > tmp
printf '@HD\tVN:1.5SO:coordinate\n' > "${sample_id}".sam
awk '/^>/{if (l!="") print l; print; l=0; next}{l+=length($0)}END{print l}' "${sample_id}"_uniref90_unaligned.fa | awk '{ print $1 }' | awk '{$1=$1;printf("%s ",$0)};NR%2==0{print ""}' | sed 's/ *$//' | sed 's/>/@SQ\tSN:/g' | sed 's/ /\tLN:/g' | awk '!x[$2]++' >> "${sample_id}".sam
head -n 3 tmp | tail -n 2 >> "${sample_id}".sam
awk 'BEGIN {FS=OFS="\t"} {$2=3; $3=$1; $4=1; $5=100; l=length($10); m="M"; s=l m ; $6=s; print}' tmp | tail -n +4 >> "${sample_id}".sam
samtools view -bS "${sample_id}".sam | samtools sort > "${sample_id}".unaligned.bam && samtools index "${sample_id}".unaligned.bam && rm "${sample_id}".sam
