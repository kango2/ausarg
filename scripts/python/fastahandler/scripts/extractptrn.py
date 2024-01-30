# FASTAhandler extran pattern match seqeucne created by Hyungtaek Jung
# A multi-fasta (multiline) file to find, filter and extract based on different criteria (filter length and seqeuence patterns).
# Example command: python3 extractptrn.py --input-seq test.fa --input-ptrn find_ptrn.txt --len-over 100 --out test_out.fa (w/ optional for --t cpu and --mem memory) 


#!/usr/bin/env python3

import argparse
import os
import re
import gzip
import logging
from collections import defaultdict

# Set up logging
logging.basicConfig(filename='extractptr.log', level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')

def unzip_file(input_file, output_dir):
    # Check if the file is compressed and extract if necessary
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
    sequences = {seq_id: sequence.replace(" ", "") for seq_id, sequence in sequences.items()}
    return sequences

def filter_by_pattern(sequences, pattern_file):
    if not pattern_file:
        return sequences

    with open(pattern_file, "r") as file:
        patterns = set(line.strip() for line in file if line.strip())

    filtered_sequences = {}
    for seq_id, sequence in sequences.items():
        for pattern in patterns:
            if pattern in sequence:
                filtered_sequences[seq_id] = sequence
                break

    return filtered_sequences

def filter_by_length(sequences, length_threshold):
    if not length_threshold:
        return sequences

    filtered_sequences = {seq_id: sequence for seq_id, sequence in sequences.items() if len(sequence) >= length_threshold}
    return filtered_sequences

def generate_fasta_output(sequences, output_file):
    with open(output_file, 'w') as file:
        for seq_id, sequence in sequences.items():
            file.write(f">{seq_id}\n{sequence}\n")

    logging.info(f"Filtered fasta output written to {output_file} successfully.")

def main():
    parser = argparse.ArgumentParser(description="Comprehensive Python script for extracting patterns from fasta files.")
    parser.add_argument("--input-seq", required=True, help="Input fasta file path.")
    parser.add_argument("--input-ptrn", help="Input pattern file path.")
    parser.add_argument("--len-over", type=int, help="Length to filter out sequences.")
    parser.add_argument("--out", required=True, help="Output fasta file path.")
    parser.add_argument("--t", type=int, default=1, help="Number of CPUs.")
    parser.add_argument("--mem", type=int, default=10, help="Amount of memory in gigabytes.")

    args = parser.parse_args()

    # Check for valid CPU and memory values
    if args.t < 1 or args.mem < 1:
        raise ValueError("Not an integer number. Please provide only readable integer numbers to call the proper CPUs and Memory.")
    logging.info(f"Using {args.t} CPUs and {args.mem}GB of memory.")

    # Check and unzip input sequence file if necessary
    if not os.path.isfile(args.input_seq):
        raise FileNotFoundError("Input sequence file not found.")
    logging.info(f"Input sequence file: {args.input_seq}")

    # Parse fasta sequences from input file
    sequences = parse_fasta(args.input_seq)

    # Filter sequences based on input pattern file
    sequences = filter_by_pattern(sequences, args.input_ptrn)

    # Filter sequences based on input length threshold
    sequences = filter_by_length(sequences, args.len_over)

    # Generate filtered fasta output
    generate_fasta_output(sequences, args.out)

if __name__ == "__main__":
    main()

