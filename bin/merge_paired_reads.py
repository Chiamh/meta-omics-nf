import sys
import gzip
from Bio import SeqIO

def re_pair_fastq_reads(pan_r1_file, pan_r2_file, dmnd_r1_file, dmnd_r2_file, sample_id):
    """
    Re-pairs FASTQ reads from two files based on read IDs and writes paired reads
    to new files and singletons to a separate file.
    """
    output_paired_r1,output_paired_r2=sample_id+"_pan_dmnd_merged_1.fastq",sample_id+"_pan_dmnd_merged_2.fastq"
    paired_r1,paired_r2 = [],[]

    for r1_file,r2_file in [(pan_r1_file,pan_r2_file),(dmnd_r1_file,dmnd_r2_file)]:
        r1_records = SeqIO.to_dict(SeqIO.parse(gzip.open(r1_file,'rt'), "fastq"))
        r2_records = SeqIO.to_dict(SeqIO.parse(gzip.open(r2_file,'rt'), "fastq"))

        # Iterate through R1 records and find corresponding R2 records
        for r1_id, r1_rec in r1_records.items():
            # Clean up read IDs by removing /1 or /2 suffixes if present
            cleaned_r1_id = r1_id.split('/')[0]
            if cleaned_r1_id in r2_records:
                paired_r1.append(r1_rec)
                paired_r2.append(r2_records[cleaned_r1_id])
                del r2_records[cleaned_r1_id] # Remove from R2 dict to avoid re-matching

        # Add remaining R2 records (which are singletons)
        singletons.extend(r2_records.values())

        # Write out the results
        SeqIO.write(paired_r1, output_paired_r1, "fastq")
        SeqIO.write(paired_r2, output_paired_r2, "fastq")

pan_aligned_fastq_1,pan_aligned_fastq_2=sys.argv[1],sys.argv[2]
dmnd_aligned_fastq_1,dmnd_aligned_fastq_2=sys.argv[3],sys.argv[4]
sample_id=sys.argv[5]

re_pair_fastq_reads(pan_aligned_fastq_1, pan_aligned_fastq_2, dmnd_aligned_fastq_1, dmnd_aligned_fastq_2, sample_id)
