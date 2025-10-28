#!/bin/bash
set -e
set -u

sample_id=${1}
pangenome_annots=${2}

samtools view -F 4 "${sample_id}"_bt2_pangenome_aligned.bam | cut -f 1,3 | sort -k1,1 > "${sample_id}"_read_to_pangene.tsv

grep ">" "${sample_id}"_uniref90_unaligned.fa | sed 's/>//g' | cut -d " " -f1 | sort > "${sample_id}"_unannotated_readnames

cut -f2,3 "${sample_id}"_kraken2.out | sed 's/(taxid/taxid/g' | awk -F "taxid" -v OFS="\t" '{print $1}' | sort -k1,1 > "${sample_id}"_k2_read_sorted

join -1 1 -2 1 -e '-' -o auto -t $'\t' "${sample_id}"_k2_read_sorted "${sample_id}"_read_to_pangene.tsv | sort -t $'\t' -k3,3 \
| join -a 1 -1 3 -2 1 -t $'\t' -e '-' -o auto - <(awk -F "\t" -v OFS="\t" '(NR>1){print $1,$3}' "${pangenome_annots}") > "${sample_id}"_pangene_taxonomy.tsv

join -1 1 -2 1 -e '-' -o auto -t $'\t' "${sample_id}"_k2_read_sorted <(cut -f1,2 "${sample_id}"_uniref90_aligned.out | sort -k1,1) \
| awk -F "\t" -v OFS="\t" '{$1="-"FS$1}1' > "${sample_id}"_transl_search_taxonomy.tsv

cat "${sample_id}"_pangene_taxonomy.tsv "${sample_id}"_transl_search_taxonomy.tsv > "${sample_id}"_all_aligned_taxonomy.tsv

sed -i "s/ \t/\t/" "${sample_id}"_all_aligned_taxonomy.tsv

cut -f 1,3,4 "${sample_id}"_all_aligned_taxonomy.tsv | sed 's/ /_/g' |sort | uniq -c | sed 's/^ *//g' | sed 's/ /\t/g' > "${sample_id}"_all_aligned_taxonomy_summary.tsv

join -1 1 -2 1 -e '-' -o auto -t $'\t' "${sample_id}"_k2_read_sorted "${sample_id}"_unannotated_readnames > "${sample_id}"_unaligned_taxonomy.tsv

#remove temp files
rm "${sample_id}"_read_to_pangene.tsv
rm "${sample_id}"_k2_read_sorted
rm "${sample_id}"_pangene_taxonomy.tsv
rm "${sample_id}"_transl_search_taxonomy.tsv
rm "${sample_id}"_unannotated_readnames
