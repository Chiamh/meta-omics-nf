#!/bin/bash
set -e
set -u

sample_id=${1}

#Create some temporary files from the output of emapper.py
grep -v "^#" "${sample_id}".emapper.annotations | awk -F "\t" -v OFS="\t" '{print $6}' | sed 's,/,_,g' > "${sample_id}"_eggnog_index
grep -v "^#" "${sample_id}".emapper.annotations | awk -F "\t" -v OFS="\t" '{print $5}' | sed 's,/,_,g' > "${sample_id}"_eggnog_OGs
grep -v "^#" "${sample_id}".emapper.annotations | awk -F "\t" -v OFS="\t" '{print $1,$5,$6,$7,$10}' > "${sample_id}"_uniref90IDs