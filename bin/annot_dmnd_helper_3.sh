#!/bin/bash
set -e
set -u

sample_id=${1}
eggnog_OG_annots=${2}
uniref90_GO=${3}

#Join the results together for an annotation file linking each Uniref90 ID to the Eggnog OG info
#The *_max_annot_OGs file is the output from the extract_OG.sh script
paste "${sample_id}"_uniref90IDs "${sample_id}"_max_annot_OGs \
| sort -t $'\t' -k 6,6 | join -a 1 -1 6 -2 1 -t $'\t' - "${eggnog_OG_annots}" | sort -t $'\t' -k 2,2 > "${sample_id}"_uniref90_OGs
			
#from the translated search result file, extract the corresponding uniref90 protein AA lengths (only the hits)

awk -v OFS="\t" '{print $2,$12}' "${sample_id}"_uniref90_aligned.out | sort -k 1,1 | uniq > "${sample_id}".sizes

#Convert the dmnd result (1 based coordinates) into bed format (0 based coordinates)
#To do so, extract columns for sseqid (2nd column),  sstart -1 (10th column values -1) and ssend (11th column)
	
awk -v OFS="\t" '{print $2,$10-1,$11}' "${sample_id}"_uniref90_aligned.out | sort -k 1,1 > "${sample_id}".bed

#Run bedtools genomecov with the bed file, to obtain subject (uniref90 cluster) read coverage
# -max in bedtools genomecov combines all positions with a depth >= max into a single bin in the histogram

bedtools genomecov -i "${sample_id}".bed -g "${sample_id}".sizes -max 1 | grep -v "genome" > "${sample_id}"_subjectcov.tsv

#As a separate file, keep only the uniref90 clusters that are >= 50% covered by "well aligned reads", as proposed by Beghini et al 2021 PMID: 33944776
#Because of the -max 1 option in genomecov, here we are checking the fraction of bases on feature with depth >= to 1

awk -v OFS="\t" '($2==1 && $5 >= 0.5){$5=$5*100; print $0}' "${sample_id}"_subjectcov.tsv > "${sample_id}"_keptclusters.tsv

#Create a file for the final output and add headers. 12 columns in final output

touch "${sample_id}"_transl-search_annot.tsv

echo -e "uniref90_ID\tpercent_cov\tAA_length\tunpaired_read_count\tuniref90_desc\temapper_max_annot_OG\temapper_OG\temapper_max_annot_lvl\teggnog_cat\teggnog_GO\teggnog_desc\tuniref90_GO" >> "${sample_id}"_transl-search_annot.tsv

#Perform a series of sequential left joins to get count information and OG annotations, only for the "kept clusters"
	
awk -v OFS="\t" '{print $1,$5}' "${sample_id}"_keptclusters.tsv \
| join -a 1 -1 1 -2 1 -t $'\t' - "${sample_id}".sizes \
| join -a 1 -1 1 -2 2 -t $'\t' - <(awk -v OFS="\t" '{print $2}' "${sample_id}"_uniref90_aligned.out \
| sort | uniq -c | sed 's/^ *//; s/ /\t/') | join -a 1 -1 1 -2 1 -t $'\t' - <(cut -f 2,15 "${sample_id}"_uniref90_aligned.out | sort -t $'\t' -k 1,1 | uniq) \
| join -a 1 -1 1 -2 2 -t $'\t' -e '-' -o auto - "${sample_id}"_uniref90_OGs \
| join -a 1 -1 1 -2 1 -t $'\t' -e '-' -o auto - "${uniref90_GO}" >> "${sample_id}"_transl-search_annot.tsv