#!/bin/bash
set -e
set -u

sample_id=${1}
panalign_results=${2}
pangenome_annots=${3}

#Transfer information from pangenome annotation file to panalign results for the sample

sed '1d' "${panalign_results}" | join -a 1 -1 1 -2 1 -t $'\t' - <(sed '1d' "${pangenome_annots}") \
| sed '1 i\pangene\tlength\tpercent_cov\tunpaired_read_count\tpangene_desc\tuniref90_ID\tuniref90_desc\tuniref90_GO\temapper_max_annot_OG\temapper_OG\temapper_max_annot_lvl\teggnog_cat\teggnog_GO\teggnog_desc' > "${sample_id}"_panalign_annot.tsv