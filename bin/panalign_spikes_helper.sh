#!/bin/bash
set -e
set -u

sample_id=${1}
spike_in_path=${2}

cut -f2,3 "${sample_id}"_kraken2.out | grep -Ef ${spike_in_path} - | cut -f1 > "${sample_id}"_spike_reads

filterbyname.sh in="${sample_id}"_bt2_pangenome_aligned_unfiltered.bam out="${sample_id}"_bt2_pangenome_aligned.bam names="${sample_id}"_spike_reads include=false

filterbyname.sh in="${sample_id}"_bt2_pangenome_unaligned_unfiltered.fastq.gz out="${sample_id}"_bt2_pangenome_unaligned.fastq.gz names="${sample_id}"_spike_reads include=false

pileup.sh in="${sample_id}"_bt2_pangenome_aligned.bam out=stdout.tsv | awk -F "\t" -v OFS="\t" '(NR >1 && $5 >= 50){print $1,$3,$5,$7+$8}' \
| sort -t $'\t' -k 1,1 | sed '1 i\pangene\tlength\tpercent_cov\tunpaired_read_count' > "${sample_id}"_bt2_pangenome_aligned_filtered_cov.tsv