#!/bin/bash
set -e
set -u

sample_id=${1}
fastq_1=${2}
fastq_2=${3}
human_pangenome_path=${4}
cpus=${5}



vg giraffe -p -t "${cpus}" -Z "${human_pangenome_path}"/hprc-v1.1-mc-chm13.gbz -d "${human_pangenome_path}"/hprc-v1.1-mc-chm13.dist -m "${human_pangenome_path}"/hprc-v1.1-mc-chm13.shortread.withzip.min -z "${human_pangenome_path}"/hprc-v1.1-mc-chm13.shortread.zipcodes -f "${fastq_1}" -f "${fastq_2}" > "${sample_id}"_hprc-v1.1-mc-chm13.gam
vg stats -a "${sample_id}"_hprc-v1.1-mc-chm13.gam > "${sample_id}"_hprc-v1.1-mc-chm13_QC.txt
vg filter -U -P "${sample_id}"_hprc-v1.1-mc-chm13.gam | vg view -X - | awk 'NR%4==1 {print substr($1,2)}' > read_names
seqtk subseq "${fastq_1}" read_names > "${sample_id}"_pancleaned_1.fastq && gzip "${sample_id}"_pancleaned_1.fastq
seqtk subseq "${fastq_2}" read_names > "${sample_id}"_pancleaned_2.fastq && gzip "${sample_id}"_pancleaned_2.fastq
