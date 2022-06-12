#!/bin/bash
set -e
set -u

sample_id=${1}
unaligned_fa=${2}
k2_out=${3}
pangenome_annots=${4}
dmnd_aligned_results=${5}

grep ">" "${unaligned_fa}" | sed 's/>//g' | cut -d " " -f1 | sort > "${sample_id}"_unannotated_readnames

cut -f2,3 "${k2_out}" | sed 's/(taxid/taxid/g' | awk -F "taxid" -v OFS="\t" '{print $1}' | sort -k1,1 \
| join -a 2 -1 1 -2 1 -e '-' -o auto -t $'\t' - "${sample_id}"_read_to_pangene.tsv | sort -t $'\t' -k3,3 \
| join -a 1 -1 3 -2 1 -t $'\t' -e '-' -o auto - <(awk -F "\t" -v OFS="\t" '(NR>1){print $1,$3}' "${pangenome_annots}") > "${sample_id}"_pangene_taxonomy.tsv

cut -f2,3 "${k2_out}" | sed 's/(taxid/taxid/g' | awk -F "taxid" -v OFS="\t" '{print $1}' | sort -k1,1 \
| join -a 2 -1 1 -2 1 -e '-' -o auto -t $'\t' - <(cut -f1,2 "${dmnd_aligned_results}" | sort -k1,1) \
| awk -F "\t" -v OFS="\t" '{$1="-"FS$1}1' > "${sample_id}"_transl_search_taxonomy.tsv

cat "${sample_id}"_pangene_taxonomy.tsv "${sample_id}"_transl_search_taxonomy.tsv > "${sample_id}"_all_aligned_taxonomy.tsv

cut -f2,3 "${k2_out}" | sed 's/(taxid/taxid/g' | awk -F "taxid" -v OFS="\t" '{print $1}' | sort -k1,1 \
| join -a 2 -1 1 -2 1 -e '-' -o auto -t $'\t' - "${sample_id}"_unannotated_readnames > "${sample_id}"_unaligned_taxonomy.tsv