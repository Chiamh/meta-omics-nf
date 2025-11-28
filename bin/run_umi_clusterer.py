#!/usr/bin/env python3
"""
run_umi_clusterer.py

Process a tab-separated file with columns:
    barcode    cluster_id    count

and cluster UMIs per cluster_id using UMI-tools' UMIClusterer.

Usage:
    python run_umi_clusterer.py <input.tsv> <output.tsv> <method>

Example:
    python run_umi_clusterer.py barcodes.tsv clustered.tsv directional
"""

import sys
import csv
from collections import defaultdict
from umi_tools.network import UMIClusterer


def run_umi_clusterer(input_tsv, output_tsv, method="directional", threshold=1):
    # Group barcodes by cluster ID
    clusters = defaultdict(dict)
    with open(input_tsv, newline='') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if len(row) < 3:
                continue  # skip malformed rows
            barcode, cluster_id, count = row[0], row[1], int(row[2])
            clusters[cluster_id][barcode] = count

    clusterer = UMIClusterer(cluster_method=method)

    with open(output_tsv, 'w', newline='') as out:
        writer = csv.writer(out, delimiter='\t')
        writer.writerow(["cluster_id", "clustered_group", "UMIs_in_group"])

        for cluster_id, barcode_counts in clusters.items():
            # Encode keys to bytes for umi_tools
            barcode_counts_bytes = {k.encode(): v for k, v in barcode_counts.items()}
            #grouped is a list of lists e.g. [[b'ATAT', b'GTAT'], [b'CCAT']]
            #if you were in the business of deduplicating reads, you’d keep one read associated with the ATAT and CCAT UMIs, and discard the reads associated with GTAT
            grouped = clusterer(barcode_counts_bytes, threshold=threshold)

            # Decode groups and write only the first UMI per group (de-duplication)
            for group in grouped:
                decoded = [umi.decode() for umi in group]
                writer.writerow([cluster_id, decoded[0], len(decoded)])
                #writer.writerow([cluster_id, ",".join(decoded), len(decoded)]) only if you want all the barcodes returned as a comma separated vector in the output


if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: run_umi_clusterer.py <input.tsv> <output.tsv> <method>")
        sys.exit(1)

    input_tsv = sys.argv[1]
    output_tsv = sys.argv[2]
    method = sys.argv[3]

    run_umi_clusterer(input_tsv, output_tsv, method)
