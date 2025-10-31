#!/bin/bash

#Script to extract and summarize MTX(METATRANSCRIPTOMIC) read count information from the nextflow MTX pipeline
#Loop through a file of Library IDs e.g. MHS300
#Run this from the scripts folder

#Create output file if it does not already exist

LIBID=$1
DATADIR=$2
OUTDIR=$3

if [[ ! -e "$OUTDIR"/UMI_MTX_QC_stats.txt ]]; then
    touch "$OUTDIR"/UMI_MTX_QC_stats.txt
	echo -e "LIBID\tBEFORE_FASTP\tAFTER_FASTP\tDUP_RATE\tAFT_HUMAN_RM\tAFT_RIBO_RM\tBT2_READS_START\tAFT_DEDUP\tPANGENE_ALIGN_READ_PAIRS\tUNIREF_ALIGN_READ_PAIRS\tOVERALL_ANNOT_RATE\tK2_UNCLASSIFIED\tK2_CLASSIFIED\tMICROBE_COUNT\tBACTERIA\tARCHAEA\tFUNGI\tVIRUS" > "$OUTDIR"/UMI_MTX_QC_stats.txt
fi

echo "Processing $1"

#Stats from fastp and bbduk are not number of read pairs, but are the sum of R1 and R2. Therefore, the formated values are divided by two for read paits
#-z returns true if a bash variable is unset
#Counts of read pairs before decont, but after fastp QC
BEFORE_FASTP=`awk -F ":" '/"before_filtering"/{ getline; print $2 }' $DATADIR/decont/RNA/"$LIBID"*.json | sed 's/,//g'`

if [[ ! -z "$BEFORE_FASTP" ]]; then
BEFORE_FASTP_F=`echo $BEFORE_FASTP / 2| bc`
else
#empty value
BEFORE_FASTP_F=
fi

#Counts of read pairs after fastp QC
AFT_FASTP=`awk -F ":" '/"after_filtering"/{ getline; print $2 }' $DATADIR/decont/RNA/"$LIBID"*.json | sed 's/,//g'`

if [[ ! -z "$AFT_FASTP" ]]; then
AFT_FASTP_F=`echo $AFT_FASTP / 2| bc`
else
#empty value
AFT_FASTP_F=
fi

#Duplication rates reported by fastp
DUP=`awk -F ": " '/duplication/{ getline; print $2 }' $DATADIR/decont/RNA/"$LIBID"*.json`

#Counts after removal of human transcripts
AFT_HG_MAP=`grep "#Total" $DATADIR/decont/RNA/"$LIBID"*_rRNAfilter.log | cut -f 2`

if [[ ! -z "$AFT_HG_MAP" ]]; then
AFT_HG_MAP_F=`echo $AFT_HG_MAP / 2| bc`
else
#empty value
AFT_HG_MAP_F=
fi

#reads matching rRNAs
RIBORNA_COUNT=`grep "#Matched" $DATADIR/decont/RNA/"$LIBID"*_rRNAfilter.log | cut -f 2`

#Counts after removal of rRNAs
AFT_RIBORNA=`echo "$AFT_HG_MAP - $RIBORNA_COUNT" | bc`

if [[ ! -z "$AFT_RIBORNA" ]]; then
AFT_RIBORNA_F=`echo $AFT_RIBORNA / 2| bc`
else
#empty value
AFT_RIBORNA_F=
fi


#Unpaired Reads entering bowtie2 for pangene alignment, before de-duplication
BT_START=`grep "were unpaired" $DATADIR/MTX_panalign_out/"$LIBID"*_bt2.log | sed 's/^ *//g' | cut -d ' ' -f 1`

#Convert to number of paired reads
if [[ ! -z "$BT_START" ]]; then
BT_START_F=`echo $BT_START / 2| bc`
else
#empty value
BT_START_F=
fi

#Kraken2 unclassified read count (paired end)
UNCLASSIFIED=`awk '($5==0){print $2}' $DATADIR/kraken2/RNA/"$LIBID"*_kraken2.tax`

#Kraken2 classified read count (paired end), classified to Root, not necessarily at species level
CLASSIFIED=`awk '($5==1){print $2}' $DATADIR/kraken2/RNA/"$LIBID"*_kraken2.tax`

#Counts after UMI based dedup equals the total number of reads entering Kraken2 = sum of classified and unclassified read pairs
AFT_DEDUP_F=`echo "$UNCLASSIFIED + $CLASSIFIED" | bc`

if [[ -z "$AFT_DEDUP_F" ]]; then
#empty value
AFT_DEDUP_F=
fi

#Alignment rates to microbial gene and protein sequences

#Number of read pairs with at least one member of the pair aligning to microbial pangenes, after de-duplication (if enabled). Includes low coverage features.
PANGENE_ALIGN_READS=`grep "_pan_umi_read_names" $DATADIR/decont/RNA/"$LIBID"*_all_dedup_read_pair_stats.txt | awk '{print $1}'`

#Number of read pairs with at least one member of the pair aligning to Uniref90, after de-duplication (if enabled). Includes low coverage features.
UNIREF_ALIGN_READS=`grep "_dmnd_umi_read_names" $DATADIR/decont/RNA/"$LIBID"*_all_dedup_read_pair_stats.txt | awk '{print $1}'`

OVERALL_ALIGN_READS=`echo $PANGENE_ALIGN_READS + $UNIREF_ALIGN_READS| bc`

#Overall annotation rate to pangene and uniref90, amongst classified reads
OVERALL_ALIGN_RATE=`echo "scale=4; ($OVERALL_ALIGN_READS / $CLASSIFIED) * 100" | bc`


#Left over reads classified as "Homo sapiens" by kraken2 
#HUMAN=`grep "Homo sapiens$" $DATADIR/kraken2/RNA/"$LIBID"_merged_kraken2.tax | cut -f2`

#Fungal reads, not necessarily at species level
FUNGI=`awk '($5==4751){print $2}' $DATADIR/kraken2/RNA/"$LIBID"*_kraken2.tax`

#Bacteria reads, not necessarily at species level
BACTERIA=`awk '($5==2){print $2}' $DATADIR/kraken2/RNA/"$LIBID"*_kraken2.tax`

#Archaea reads, not necessarily at species level
ARCHAEA=`awk '($5==2157){print $2}' $DATADIR/kraken2/RNA/"$LIBID"*_kraken2.tax`

#Viral reads, not necessarily at species level
VIRUS=`awk '($5==10239){print $2}' $DATADIR/kraken2/RNA/"$LIBID"*_kraken2.tax`

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


#Microbial PAIRED read count (after human, rRNA removal and de-duplication). Excludes "unclassified"
MICROBE_COUNT=`echo "$FUNGI + $BACTERIA + $ARCHAEA + $VIRUS" | bc`

#Outputs
#echo is needed for printing in new line
echo -e "$LIBID\t$BEFORE_FASTP_F\t$AFT_FASTP_F\t$DUP\t$AFT_HG_MAP_F\t$AFT_RIBORNA_F\t$BT_START_F\t$AFT_DEDUP_F\t$PANGENE_ALIGN_READS\t$UNIREF_ALIGN_READS\t$OVERALL_ALIGN_RATE\t$UNCLASSIFIED\t$CLASSIFIED\t$MICROBE_COUNT\t$BACTERIA\t$ARCHAEA\t$FUNGI\t$VIRUS" >> "$OUTDIR"/UMI_MTX_QC_stats.txt
