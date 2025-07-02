#!/bin/bash
set -e
set -u

bam=${1}

pileup.sh in="${bam}" out=stdout.tsv | awk -F "\t" -v OFS="\t" '(NR >1 && $5 >= 50){print $1,$3,$5,$7+$8}' \
| sort -t $'\t' -k 1,1 | sed '1 i\pangene\tlength\tpercent_cov\tunpaired_read_count' > "${bam}"_cov.tsv
