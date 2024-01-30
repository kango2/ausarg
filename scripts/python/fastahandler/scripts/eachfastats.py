# FASTAhandler each ID and length created by Hyungtaek Jung
# A multi-fasta (multiline) file to calculate each fasta sequence stats
# Example command: python3 eachfastats.py --input-seq test1.fa --out test_stats_out.txt (w/ optional for --t cpu and --mem memory)


#!/usr/bin/env python3

import argparse
import os
import re
import gzip
import logging
from collections import defaultdict

# Set up logging
logging.basicConfig(filename='eachfastats.log', level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')

def unzip_file(input_file, output_dir):
    # Check compressed file and extract it
    if input_file.lower().endswith(('.tar.gz', '.gz', '.zip', '.bz2')):
        with gzip.open(input_file, 'rb') as f_in:
            with open(os.path.join(output_dir, os.path.basename(input_file)[:-3]), 'wb') as f_out:
                f_out.write(f_in.read())
        return os.path.join(output_dir, os.path.basename(input_file)[:-3])
    return input_file

def parse_fasta(file_path):
    # Check compressed file and extract it
    file_path = unzip_file(file_path, os.path.dirname(file_path))

    sequences = defaultdict(str)
    current_seq_id = None
    with open(file_path, "r") as file:
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                current_seq_id = line[1:]
            else:
                sequences[current_seq_id] += line

    # Convert multiline sequences to single-line format
    sequences = {seq_id: sequence.replace(" ", "").upper() for seq_id, sequence in sequences.items()}
    return sequences

def calculate_summary_stats(file_path):
    sequences = parse_fasta(file_path)
    summary_stats = [(re.sub(r'>', '', seq_id), len(re.sub(r'[^ATGCN]', '', sequence))) for seq_id, sequence in sequences.items()]
    return summary_stats

def generate_summary_output(summary_stats, output_file):
    with open(output_file, 'w') as file:
        for input_id, sequence_length in summary_stats:
            file.write(f"{input_id}\t{sequence_length}\n")

    logging.info(f"Summary output written to {output_file} successfully.")

def main():
    parser = argparse.ArgumentParser(description="Comprehensive Python script for analyzing fasta files.")
    parser.add_argument("--input-seq", required=True, help="Input fasta file path.")
    parser.add_argument("--out", required=True, help="Output file path for summary results.")
    parser.add_argument("--t", type=int, default=1, help="Number of CPUs.")
    parser.add_argument("--mem", type=int, default=10, help="Amount of memory in gigabytes.")

    args = parser.parse_args()

    # Check for valid CPU and memory values
    if args.t < 1 or args.mem < 1:
        raise ValueError("Not an integer number. Please provide only readable integer numbers to call the proper CPUs and Memory.")
    logging.info(f"Using {args.t} CPUs and {args.mem}GB of memory.")

    # Check and unzip input file if necessary
    if not os.path.isfile(args.input_seq):
        raise FileNotFoundError("Input file not found.")
    logging.info(f"Input file: {args.input_seq}")

    # Calculate summary statistics
    summary_stats = calculate_summary_stats(args.input_seq)

    # Generate log file
    logging.info(f"Summary statistics: {summary_stats}")

    # Generate summary output
    generate_summary_output(summary_stats, args.out)

if __name__ == "__main__":
    main()

