#!/bin/bash
set -e
set -u

sample_id=${1}

pileup.sh in="${sample_id}"_bt2_pangenome_aligned.bam out=stdout.tsv | awk -F "\t" -v OFS="\t" '(NR >1 && $5 >= 50){print $1,$3,$5,$7+$8}' \
| sort -t $'\t' -k 1,1 | sed '1 i\pangene\tlength\tpercent_cov\tunpaired_read_count' > "${sample_id}"_bt2_pangenome_aligned_filtered_cov.tsv
