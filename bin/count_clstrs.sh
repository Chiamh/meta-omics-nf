#!/bin/bash

set -e
set -u

INPUT=${1}

sed 's/cluster //g' ${INPUT} | awk -v OFS="\t" '{
	key = $1 "\t" $2  #Use barcode + cluster as a key
	count[key]++
}
END { for (k in count)
	print k, count[k]
}' | sort -k2,2n -k3,3nr   #the second column is the cluster number, third column is the number of times a barcode appears in that cluster
