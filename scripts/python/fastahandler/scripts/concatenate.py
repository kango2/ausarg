# FastaHandler concatenate ID and sequence created by Hyungtaek Jung
# Concatenate sequences with the same ID from multiple fasta files (unlimited)
# Example command: python3 concatenate.py --input-seq test1.fa test2.fa test3.fa test4.fa test5.fa --out concat_seq_out.fa (w/ optional for --t cpu and --mem memory)


#!/usr/bin/env python3

import argparse
import os
import gzip
import logging
from collections import defaultdict

# Set up logging
logging.basicConfig(filename='concatenate.log', level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')

def unzip_file(input_file, output_dir):
    # Check if the file is compressed and extract if necessary
    compression_extensions = ('.tar.gz', '.gz', '.zip', '.bz2')
    if input_file.lower().endswith(compression_extensions):
        with gzip.open(input_file, 'rb') as f_in:
            with open(os.path.join(output_dir, os.path.basename(input_file)[:-3]), 'wb') as f_out:
                f_out.write(f_in.read())
        return os.path.join(output_dir, os.path.basename(input_file)[:-3])
    return input_file

def parse_fasta(file_path):
    # Check compressed and extract it
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

    return sequences

def concatenate_sequences(input_files):
    if not input_files:
        raise ValueError("At least one input sequence file is required.")

    all_sequences = defaultdict(str)
    for file in input_files:
        sequences = parse_fasta(file)
        for seq_id, sequence in sequences.items():
            all_sequences[seq_id] += sequence

    return all_sequences

def generate_fasta_output(sequences, output_file):
    with open(output_file, 'w') as file:
        for seq_id, sequence in sequences.items():
            file.write(f">{seq_id}\n{sequence}\n")

    logging.info(f"Concatenated fasta output written to {output_file} successfully.")

def main():
    parser = argparse.ArgumentParser(description="Comprehensive Python script for concatenating multiple fasta files.")
    parser.add_argument("--input-seqs", nargs='+', required=True, help="Input fasta file paths for sequences.")
    parser.add_argument("--out", required=True, help="Output fasta file path.")
    parser.add_argument("--t", type=int, default=1, help="Number of CPUs.")
    parser.add_argument("--mem", type=int, default=10, help="Amount of memory in gigabytes.")

    args = parser.parse_args()

    # Check for valid CPU and memory values
    args.t = max(1, args.t)
    args.mem = max(1, args.mem)
    logging.info(f"Using {args.t} CPUs and {args.mem}GB of memory.")

    # Check and unzip input sequence files if necessary
    input_files = args.input_seq
    for file in input_files:
        if not os.path.isfile(file):
            raise FileNotFoundError(f"Input sequence file not found: {file}")
        logging.info(f"Input sequence file: {file}")

    # Concatenate sequences from multiple input files
    concatenated_sequences = concatenate_sequences(input_files)

    # Generate concatenated fasta output
    generate_fasta_output(concatenated_sequences, args.out)

if __name__ == "__main__":
    main()

