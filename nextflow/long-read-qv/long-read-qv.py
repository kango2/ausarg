import sys
import os
from Bio import SeqIO
import csv
import gzip

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

def process_fastq(input_fastq, output_path):
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

    input_basename = os.path.splitext(os.path.basename(input_fastq))[0]
    quality_output_csv = os.path.join(output_path, f"{input_basename}_quality_freq.csv")
    length_output_csv = os.path.join(output_path, f"{input_basename}_length_freq.csv")
    stats_output_csv = os.path.join(output_path, f"{input_basename}_stats.csv")

    # Write quality frequency data
    with open(quality_output_csv, "w", newline="") as csvfile:
            csv_writer = csv.writer(csvfile)
            csv_writer.writerow(["File_Path", "Read_Length", "QV", "Read_Numbers"])
            for bin_key, bin_data in bins.items():
                length_bin, qv_bin = bin_key
                frequency = bin_data["count"]
                csv_writer.writerow([input_fastq, length_bin, qv_bin, frequency])

    # Write length frequency data
    with open(length_output_csv, "w", newline="") as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(["File_Path", "Read_Length", "Summed_Read_Numbers"])
        for length_bin, count in length_sums.items():
            csv_writer.writerow([input_fastq, length_bin, count])


    # Write stats data
    with open(stats_output_csv, "w", newline="") as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(["File_Path", "Total_Bases", "Total_Reads", "Average_Read_Length", "N50", "N90", "L50", "L90", "Total_Ns"])
        csv_writer.writerow([input_fastq, total_bases, total_reads, average_read_length, n50, n90, l50, l90, total_ns])

    print("CSV output written to:", quality_output_csv)
    print("CSV output written to:", length_output_csv)
    print("CSV output written to:", stats_output_csv)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script_name.py input_fastq.gz output_path")
        sys.exit(1)

    input_fastq = sys.argv[1]
    output_path = sys.argv[2]

    process_fastq(input_fastq, output_path)
