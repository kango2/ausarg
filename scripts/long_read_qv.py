#!/usr/bin/env python
from Bio import SeqIO
import sys
import os
import csv
import gzip
import argparse

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

def process_fastq(input_fastq, sample, output_folder):
    bins = {}
    length_sums = {}
    read_lengths = []
    total_bases = 0
    total_reads = 0
    total_ns = 0

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

    # Ensure output directory exists
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Define file paths
    file_prefix = os.path.join(output_folder, sample)
    quality_output_csv = f"{file_prefix}_quality_freq.csv"
    length_output_csv = f"{file_prefix}_length_freq.csv"
    stats_output_csv = f"{file_prefix}_stats.csv"

    # Write quality frequency data
    with open(quality_output_csv, "w", newline="") as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(["Sample", "Read_Length", "QV", "Read_Numbers"])
        for bin_key, bin_data in bins.items():
            length_bin, qv_bin = bin_key
            frequency = bin_data["count"]
            csv_writer.writerow([sample, length_bin, qv_bin, frequency])

    # Write length frequency data
    with open(length_output_csv, "w", newline="") as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(["Sample", "Read_Length", "Summed_Read_Numbers"])
        for length_bin, count in length_sums.items():
            csv_writer.writerow([sample, length_bin, count])

    # Write stats data
    with open(stats_output_csv, "w", newline="") as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(["Sample", "Total_Bases", "Total_Reads", "Average_Read_Length", "N50", "N90", "L50", "L90", "Total_Ns"])
        csv_writer.writerow([sample, total_bases, total_reads, average_read_length, n50, n90, l50, l90, total_ns])


def main():
    parser = argparse.ArgumentParser(description="Process a FASTQ file and output statistics and frequencies.")
    parser.add_argument("-input", help="Path to the input FASTQ file", required=True)
    parser.add_argument("-sample", help="Sample name to use as file prefix", required=True)
    parser.add_argument("-output", help="Folder to store the output CSV files", required=True)
    
    args = parser.parse_args()
    
    process_fastq(args.input, args.sample, args.output)

if __name__ == "__main__":
    main()
