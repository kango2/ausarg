import sys
import os
from Bio import SeqIO
import gzip
import csv

bin_size = 100

def process_fastq(input_fastq, output_file_name):
    bins = {}

    with open(input_fastq, "rt") as f:
        for record in SeqIO.parse(f, "fastq-sanger"):
            sequence_length = len(record.seq)
            avg_qv = round(sum(record.letter_annotations["phred_quality"]) / sequence_length)

            bin_number = (sequence_length) // bin_size * bin_size

            bin_key = (bin_number, avg_qv)

            if bin_key not in bins:
                bins[bin_key] = {"total_qv": 0, "count": 0}
            bins[bin_key]["total_qv"] += avg_qv
            bins[bin_key]["count"] += 1

    output_csv = f"{output_file_name}_{os.path.splitext(os.path.basename(input_fastq))[0]}_quality_freq.csv"

    with open(output_csv, "w", newline="") as csvfile:
        csv_writer = csv.writer(csvfile)

        csv_writer.writerow(["Read_Length", "QV", "Read_Numbers"])

        for bin_key, bin_data in bins.items():
            length_bin, qv_bin = bin_key
            frequency = bin_data["count"]
            csv_writer.writerow([length_bin, qv_bin, frequency])

    print("CSV output written to:", output_csv)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script_name.py input_fastq.gz output_file_name")
        sys.exit(1)

    input_fastq = sys.argv[1]
    output_file_name = sys.argv[2]

    process_fastq(input_fastq, output_file_name)
