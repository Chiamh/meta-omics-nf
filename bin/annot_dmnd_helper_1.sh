#!/bin/bash
set -e
set -u

sample_id=${1}
uniref90_fasta=${2}

#Extract uniref90 IDs (column 15) for OG labelling by eggnog and fix formatting.

awk -F "\t" '{print $15}' "${sample_id}"_uniref90_aligned.out | sort | uniq > "${sample_id}"_clust_to_pull.txt
filterbyname.sh in="${uniref90_fasta}" out="${sample_id}"_chosen.fa names="${sample_id}"_clust_to_pull.txt include=t fixjunk
