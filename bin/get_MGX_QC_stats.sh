#!/bin/bash

#Script to extract and summarize MGX (METAGENOMIC) read count information from the nextflow MTX pipeline
#Loop through a file of Library IDs e.g. MHS400
#Run this from the scripts folder
#for i in `cat ../metadata/MGX_LIBID.txt`; do bash get_MGX_QC_stats.sh "$i" ~/local_project/skin_mtx/metadata/ || continue ; done
#Create output file if it does not already exist

LIBID=$1
DATADIR=$2
OUTDIR=$3


if [[ ! -e "$OUTDIR"/MGX_QC_stats.txt ]]; then
    touch "$OUTDIR"/MGX_QC_stats.txt
	echo -e "LIBID\tBEFORE_FASTP\tAFTER_FASTP\tDUP_RATE\tBT2_READS_START\tPANGENE_ALIGN_RATE\tPANGENE_ALIGN_READS_UNPAIRED\tUNIREF_ALIGN_READS_UNPAIRED\tK2_UNCLASSIFIED\tK2_CLASSIFIED\tMICROBE_COUNT\tBACTERIA\tARCHAEA\tFUNGI\tVIRUS" > "$OUTDIR"/MGX_QC_stats.txt
fi

echo "Processing $1"
#Outputs from fastp are not number of read pairs, but are the sum of R1 and R2. So have to divide by 2
#-z returns true if a bash variable is unset
#Counts of read pairs before decont, but after fastp QC
BEFORE_FASTP=`awk -F ":" '/"before_filtering"/{ getline; print $2 }' $DATADIR/decont/DNA/"$LIBID"*.json | sed 's/,//g'`
#Convert to number of read pairs
#https://www.cyberciti.biz/faq/unix-linux-bash-script-check-if-variable-is-empty/  (-z checks for empty variable)
if [[ ! -z "$BEFORE_FASTP" ]]; then
BEFORE_FASTP_F=`echo $BEFORE_FASTP / 2| bc`
else
#empty value
BEFORE_FASTP_F=
fi

#Counts of read pairs after fastp QC
AFT_FASTP=`awk -F ":" '/"after_filtering"/{ getline; print $2 }' $DATADIR/decont/DNA/"$LIBID"*.json | sed 's/,//g'`
if [[ ! -z "$AFT_FASTP" ]]; then
AFT_FASTP_F=`echo $AFT_FASTP / 2| bc`
else
#empty value
AFT_FASTP_F=
fi

#Duplication rates reported by fastp
DUP=`awk -F ": " '/duplication/{ getline; print $2 }' $DATADIR/decont/DNA/"$LIBID"*.json`

#Unpaired Reads entering bowtie2 for pangene alignment
#This is also the number of reads after host removal post STAR and post vg giraffe
BT_START=`grep "were unpaired" $DATADIR/MGX_panalign_out/"$LIBID"*_bt2.log | sed 's/^ *//g' | cut -d ' ' -f 1`

#Convert to number of paired reads
if [[ ! -z "$BT_START" ]]; then
BT_START_F=`echo $BT_START / 2| bc`
else
#empty value
BT_START_F=
fi

#Alignment rate to pangenes
BT_ALIGN_RATE=`grep "overall alignment rate" $DATADIR/MGX_panalign_out/"$LIBID"*_bt2.log | cut -d '%' -f 1`

#Number of unpaired reads aligning to microbial pangenes, after de-duplication (if enabled). Includes low coverage features.
PANGENE_ALIGN_SINGLEMAP=`grep "aligned exactly 1 time" $DATADIR/MGX_panalign_out/"$LIBID"*_bt2.log |tr -s ' ' |awk '{print $1}'`
PANGENE_ALIGN_MULTIMAP=`grep "aligned >1 times" $DATADIR/MGX_panalign_out/"$LIBID"*_bt2.log |tr -s ' ' |awk '{print $1}'`
PANGENE_ALIGN_READS=`echo $PANGENE_ALIGN_SINGLEMAP + $PANGENE_ALIGN_MULTIMAP| bc`

#Number of unpaired reads aligning to Uniref90, after de-duplication (if enabled). Includes low coverage features.
UNIREF_ALIGN_READS=`grep "queries aligned." $DATADIR/MGX_dmnd_out/"$LIBID"*_dmnd.log | awk '{print $1}'`

#Kraken2 unclassified read count (paired end)
UNCLASSIFIED=`awk '($5==0){print $2}' $DATADIR/kraken2/DNA/"$LIBID"*_kraken2.tax`

#Kraken2 classified read count (paired end), classified to Root, not necessarily at species level
CLASSIFIED=`awk '($5==1){print $2}' $DATADIR/kraken2/DNA/"$LIBID"*_kraken2.tax`

#Left over reads classified as "Homo sapiens" by kraken2 
#HUMAN=`grep "Homo sapiens$" $DATADIR/kraken2/DNA/"$LIBID"_merged_kraken2.tax | cut -f2`

#Fungal reads, not necessarily at species level
FUNGI=`awk '($5==4751){print $2}' $DATADIR/kraken2/DNA/"$LIBID"*_kraken2.tax`

#Bacteria reads, not necessarily at species level
BACTERIA=`awk '($5==2){print $2}' $DATADIR/kraken2/DNA/"$LIBID"*_kraken2.tax`

#Archaea reads, not necessarily at species level
ARCHAEA=`awk '($5==2157){print $2}' $DATADIR/kraken2/DNA/"$LIBID"*_kraken2.tax`

#Viral reads, not necessarily at species level
VIRUS=`awk '($5==10239){print $2}' $DATADIR/kraken2/DNA/"$LIBID"*_kraken2.tax`

if [[ -z "$FUNGI" ]]; then
FUNGI=0
fi

if [[ -z "$BACTERIA" ]]; then
BACTERIA=0
fi

if [[ -z "$ARCHAEA" ]]; then
ARCHAEA=0
fi

if [[ -z "$VIRUS" ]]; then
VIRUS=0
fi


#Microbial PAIRED read count (after human, rDNA removal and de-duplication). Excludes "unclassified"
MICROBE_COUNT=`echo "$FUNGI + $BACTERIA + $ARCHAEA + $VIRUS" | bc`


#Outputs
#echo is needed for printing in new line
echo -e "$LIBID\t$BEFORE_FASTP_F\t$AFT_FASTP_F\t$DUP\t$BT_START_F\t$BT_ALIGN_RATE\t$PANGENE_ALIGN_READS\t$UNIREF_ALIGN_READS\t$UNCLASSIFIED\t$CLASSIFIED\t$MICROBE_COUNT\t$BACTERIA\t$ARCHAEA\t$FUNGI\t$VIRUS" >> "$OUTDIR"/MGX_QC_stats.txt
