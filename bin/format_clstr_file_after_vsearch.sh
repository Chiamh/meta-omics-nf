#!/bin/bash

set -e
set -u

INPUT=${1}

awk -F "\t" -v OFS="\t" '($1 != "C"){match($9, /_([ACGTN]+)$/); barcode = substr($9, RSTART, RLENGTH); $9 = barcode; print $9, $2}' ${INPUT} | sed 's/_//g' | sort -k2,2n
