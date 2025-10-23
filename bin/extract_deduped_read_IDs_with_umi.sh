#!/bin/bash

set -e
set -u

#The BARCODE FILE has at least two columns. First column is the cluster ID (an integer value), second column is a barcode of ATGCs
BARCODE_FILE=${1}
# READ FILE has at least two columns. First column is the read ID with barcode separated by "_", second column is the cluster ID
READ_FILE=${2}
#Library ID
LIBID=${3}

#In awk, arrays are associative (or hash) arrays. That means you can use strings as keys
#map is the name of an array.
# [id, barcode] is the key. In awk, multiple values separated by a comma are concatenated internally with a special character (the ASCII SUB character, \034), so map[id, barcode] is effectively a unique key for this pair.
# = map[clust_id, barcode] stores a value (here just 1) at that key. You do not really care about the value; you just want to record that this ID+barcode pair exists.

awk -F'\t' '
# Step 1: Read file1 into a lookup array, skipping the header row
NR==FNR && FNR>1 {
    clust_id = $1
    barcode = $2
    map[clust_id, barcode] = 1
    next
}

# Step 2: Process file2
#For each row in file2.tsv:
#Extract the barcode from the read ID (everything after _)
#Check if the ID and barcode exist in file1
#Check if this barcode has already been output for this ID (!seen[id, barcode])
# If both are true, print the line and mark it as seen
#Only the first occurrence of each barcoded read ID (in column 1) per cluster ID is kept.
#you can print $0 as sanity check
{
    split($1, arr, "_")       # arr[2] is the barcode part
    barcode = arr[2]
    clust_id = $2
    if ((clust_id, barcode) in map && !seen[clust_id, barcode]) {
        print $1   
        seen[clust_id, barcode] = 1
    }
}
' "${BARCODE_FILE}" "${READ_FILE}" | sort | uniq > "${LIBID}"_read_names_from_unaligned_dedup
