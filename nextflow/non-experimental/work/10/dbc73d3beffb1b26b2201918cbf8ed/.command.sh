#!/usr/bin/env python
import subprocess
import sys
import os
from Bio import SeqIO
import csv
import gzip
import argparse
import datetime

subprocess.run(['module unload python3'])
subprocess.run(['module unload biopython'])
subprocess.run(['module load biopython'])
bin_size = 100

def calculate_n50_n90(read_lengths):
    sorted_lengths = sorted(read_lengths, reverse=True)
    total_length = sum(sorted_lengths)
    cumulative_length = 0
    n50 = n90 = l50 = l90 = 0

    for i, length in enumerate(sorted_lengths):
        cumulative_length += length
        if not n50 and cumulative_length >= total_length * 0.5:
            n50 = length
            l50 = i + 1
        if not n90 and cumulative_length >= total_length * 0.9:
            n90 = length
            l90 = i + 1
            break

    return n50, n90, l50, l90

def log_progress(message, log_file, flowcell_id, input_file):
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    with open(log_file, "a") as log:
        log.write(f"{timestamp} - {message} - Flowcell: {flowcell_id} - File: {input_file}")


def process_fastq(input_fastq, output_path, flowcell_id, platform,sample, log_file):
    bins = {}
    length_sums = {}
    read_lengths = []
    total_bases = 0
    total_reads = 0
    total_ns = 0

    log_progress("Starting FASTQ file processing.", log_file,flowcell_id, input_fastq)

    with gzip.open(input_fastq, "rt") as f:
        for record in SeqIO.parse(f, "fastq-sanger"):
            sequence_length = len(record.seq)
            avg_qv = round(sum(record.letter_annotations["phred_quality"]) / sequence_length)
            total_bases += sequence_length
            total_reads += 1
            total_ns += record.seq.count("N")
            read_lengths.append(sequence_length)

            bin_number = (sequence_length) // bin_size * bin_size
            bin_key = (bin_number, avg_qv)

            if bin_key not in bins:
                bins[bin_key] = {"total_qv": 0, "count": 0}
            bins[bin_key]["total_qv"] += avg_qv
            bins[bin_key]["count"] += 1

            if bin_number not in length_sums:
                length_sums[bin_number] = 0
            length_sums[bin_number] += 1

    n50, n90, l50, l90 = calculate_n50_n90(read_lengths)
    average_read_length = total_bases / total_reads if total_reads > 0 else 0

    log_progress("FASTQ file processing completed. Writing CSV files...", log_file, flowcell_id, input_fastq)

    file_prefix = f"{sample}_{flowcell_id}_{platform}"
    quality_output_csv = os.path.join(output_path, f"{file_prefix}_quality_freq.csv")
    length_output_csv = os.path.join(output_path, f"{file_prefix}_length_freq.csv")
    stats_output_csv = os.path.join(output_path, f"{file_prefix}_stats.csv")

    # Write quality frequency data
    with open(quality_output_csv, "w", newline="") as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(["File_Path","Sample","Flowcell_ID", "Platform", "Read_Length", "QV", "Read_Numbers"])
        for bin_key, bin_data in bins.items():
            length_bin, qv_bin = bin_key
            frequency = bin_data["count"]
            csv_writer.writerow([input_fastq,sample, flowcell_id, platform, length_bin, qv_bin, frequency])

    # Write length frequency data
    with open(length_output_csv, "w", newline="") as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(["File_Path","Sample", "Flowcell_ID", "Platform", "Read_Length", "Summed_Read_Numbers"])
        for length_bin, count in length_sums.items():
            csv_writer.writerow([input_fastq,sample, flowcell_id, platform, length_bin, count])

    # Write stats data
    with open(stats_output_csv, "w", newline="") as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(["File_Path","Sample", "Flowcell_ID", "Platform", "Total_Bases", "Total_Reads", "Average_Read_Length", "N50", "N90", "L50", "L90", "Total_Ns"])
        csv_writer.writerow([input_fastq,sample, flowcell_id, platform, total_bases, total_reads, average_read_length, n50, n90, l50, l90, total_ns])

    log_progress("CSV files successfully written.", log_file, flowcell_id, input_fastq)
    print(f"CSV output and log file written to: {output_path}")

def main():

    file_prefix = f"{BASDU}_{PAF987}_{PACBIO_SMRT}"
    log_file = os.path.join(pipeline, f"{file_prefix}_processing_log.txt")
    process_fastq(testhifi.fq.gz, pipeline, PAF987, PACBIO_SMRT, BASDU, log_file)

if __name__ == "__main__":
    main()
