#!/bin/bash
set -e
set -u

sample_id=${1}
dmndfai=${2}

grep "UniRef90" ${sample_id}_uniref90_aligned.out | awk '{print $3}' > ${sample_id}_uniref90ids.txt
head -n 1 ${sample_id}_uniref90_aligned.out > ${sample_id}.sam
grep -Ff ${sample_id}_uniref90ids.txt ${dmndfai} | awk '{print "@SQ\tSN:" $1 "\tLN:" $2}' >> ${sample_id}.sam
tail -n +2 ${sample_id}_uniref90_aligned.out >> ${sample_id}.sam
sed -i 's/PN:DIAMOND/ID:DIAMOND/' ${sample_id}.sam
sed -i '/@mm/d' ${sample_id}.sam
