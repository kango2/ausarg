import argparse
import os
import pyfastx
import datetime

def calculate_n50_n90(read_lengths):
    sorted_lengths = sorted(read_lengths, reverse=True)
    total_length = sum(sorted_lengths)
    cumulative_length = 0
    n50 = n90 = 0

    for length in sorted_lengths:
        cumulative_length += length
        if not n50 and cumulative_length >= total_length * 0.5:
            n50 = length
        if not n90 and cumulative_length >= total_length * 0.9:
            n90 = length
            break

    return n50, n90

def log_progress(message, log_file, flowcell_id, input_file):
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    with open(log_file, "a") as log:
        log.write(f"{timestamp} - {message} - Flowcell: {flowcell_id} - File: {input_file}\n")

def process_fastq(input_fastq, output_path, flowcell_id, platform, sample, log_file):
    read_lengths = []
    total_bases = 0
    total_reads = 0

    log_progress("Starting FASTQ file processing with pyfastx.", log_file, flowcell_id, input_fastq)

    fq = pyfastx.Fastq(input_fastq, build_index=False)
    for record in fq:
        sequence_length = len(record.seq)
        total_bases += sequence_length
        total_reads += 1
        read_lengths.append(sequence_length)
    fq.close()

    n50, n90 = calculate_n50_n90(read_lengths)

    stats_output_csv = os.path.join(output_path, f"{sample}_{flowcell_id}_{platform}_stats.csv")

    # Write stats data
    with open(stats_output_csv, "w") as csvfile:
        csvfile.write("Total_Reads,Total_Bases,N50,N90\n")
        csvfile.write(f"{total_reads},{total_bases},{n50},{n90}\n")

    log_progress("FASTQ file processing completed. Stats CSV file written.", log_file, flowcell_id, input_fastq)
    print(f"Stats CSV output and log file written to: {output_path}")

def main():
    parser = argparse.ArgumentParser(description="Process a FASTQ file and generate statistics using pyfastx.")
    parser.add_argument("--i", required=True, help="Input FASTQ file (gzipped)")
    parser.add_argument("--o", required=True, help="Path to output the CSV file")
    parser.add_argument("--f", required=True, help="Flowcell ID")
    parser.add_argument("--p", required=True, help="Sequencing platform")
    parser.add_argument("--s", required=True, help="Sample name")

    args = parser.parse_args()

    file_prefix = f"{args.s}_{args.f}_{args.p}"
    log_file = os.path.join(args.o, f"{file_prefix}_processing_log.txt")
    process_fastq(args.i, args.o, args.f, args.p, args.s, log_file)

if __name__ == "__main__":
    main()
