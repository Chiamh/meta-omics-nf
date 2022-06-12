#!/bin/bash
#argument 1 is the path to a file, with each line being the "max_annot_level" (column 6) from an *.emapper.annotations file.
#argument 2 is the path to another file, with each line being the "eggNOG_OGs" (column 5) from the same *.emapper.annotations file.
#argument 3 is the sample ID
#This script extracts the OG corresponding to the "max_annot_level" for a "best hit". 
#It is the widest OG for which annotations were transferred to the queries during eggnog mapping
#output will be a single column with the "best" eggnog OG for each feature
#https://stackoverflow.com/questions/31669103/extract-part-between-two-strings-in-a-variable-in-unix

set -e
set -u

mapfile -t ARRAY < $1
count=0
while read -r line;
do
	echo "$line"| sed -e "s/.*,\(.*\)${ARRAY[$count]}.*/\1/" | cut -f 1 -d "@" >> "$3"_max_annot_OGs
       	((count+=1))
done < $2
#Note the use of double quote marks in the sed command, very important for variable expansion
