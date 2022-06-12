#!/bin/bash
set -e
set -u

sample_id=${1}

awk -F "\t" -v OFS="\t" '($4=="S" || $4 =="U"){print $1,$2,$6}' "${sample_id}"_kraken2.tax | tr -s ' ' | sed 's/^ //g' > "${sample_id}"_k2.s.tsv