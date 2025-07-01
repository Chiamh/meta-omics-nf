#!/bin/bash
set -e
set -u

sample_id=${1}
dmndfai=${2}

grep "UniRef90" ${sample_id}_uniref90_aligned.out | awk '{print $3}' > ${sample_id}_uniref90ids.txt
echo "@HD	VN:1.5	SO:query" > ${sample_id}.sam
grep -Ff ${sample_id}_uniref90ids.txt ${dmndfai} | awk '{print "@SQ\tSN:" $1 "\tLN:" $2}' >> ${sample_id}.sam
echo "@PG	ID:DIAMOND	VN:2.0.12	CL:/opt/conda/bin/diamond blastx --query MHS434_bt2_pangenome_unaligned.fastq.gz --id 80 --query-cover 90 --threads 18 --max-target-seqs 1 -b6 --outfmt 101 --db /home/users/astar/gis/wijayai/scratch/genomeDB//uniref90_09Jun2021/uniref90_09jun2021.dmnd --out MHS434_uniref90_aligned.out --un MHS434_uniref90_unaligned.fa --unfmt fasta" >> ${sample_id}.sam
echo "@CO	BlastX-like alignments" >> ${sample_id}.sam
echo "@CO	Reporting AS: bitScore, ZR: rawScore, ZE: expected, ZI: percent identity, ZL: reference length, ZF: frame, ZS: query start DNA coordinate" >> ${sample_id}.sam
tail -n +2 ${sample_id}_uniref90_aligned.out >> ${sample_id}.sam
sed -i 's/PN:DIAMOND/ID:DIAMOND/' ${sample_id}.sam
sed -i '/@mm/d' ${sample_id}.sam
