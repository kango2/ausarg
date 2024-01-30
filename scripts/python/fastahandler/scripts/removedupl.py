# FastaHandler header count and remove duplication created by Hyungtaek Jung
# Remove the duplicate ID and sequence from a multi-fasta (multiline) file to be written in a single-line fasta
# Include duplicated headers and their occurrence counts
# Example command: python3 removedupl.py --input-seq test.fa --out . (w/ optional for --t cpu and --mem memory)


#!/usr/bin/env python3

import argparse
import os
import sys
import gzip
import tarfile
import zipfile
import bz2
import logging
from collections import Counter

logging.basicConfig(filename="removedupl.log", level=logging.INFO, format="%(asctime)s - %(message)s")

def read_fasta(file_path):
    sequences = {}
    current_header = None

    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                current_header = line[1:]
                # Initialize the count for each header
                if current_header not in sequences:
                    sequences[current_header] = 1
                else:
                    sequences[current_header] += 1

    return sequences

def write_header_counts(headers, output_path):
    # Check if the output_path is a directory
    if os.path.isdir(output_path):
        output_path = os.path.join(output_path, "out_header_count.txt")

    with open(output_path, 'w') as file:
        file.write("Header\tCount\n")
        for header, count in sorted(headers.items()):
            file.write(f"{header}\t{count}\n")

def parse_args():
    parser = argparse.ArgumentParser(description="Remove duplicates from a multi-line FASTA file and count occurrences of each header.")
    parser.add_argument("--input-seq", required=True, help="Path to the input multi-line FASTA file.")
    parser.add_argument("--output", default="./", help="Path to the output file or directory for header counts.")
    parser.add_argument("--t", type=int, default=1, help="Number of CPUs with only numbers (integers). Default: 1.")
    parser.add_argument("--mem", type=int, default=10, help="Number of memory in gigabytes with only numbers (integers). Default: 10.")
    args = parser.parse_args()

    if not os.path.exists(args.input_seq):
        sys.exit("Error: Input FASTA file not found or inaccessible.")

    if not isinstance(args.t, int) or args.t <= 0:
        sys.exit("Error: Number of CPUs must be a positive integer.")

    if not isinstance(args.mem, int) or args.mem <= 0:
        sys.exit("Error: Memory must be a positive integer.")

    return args

def unzip_file(file_path):
    file_ext = os.path.splitext(file_path)[-1]
    unzip_path = file_path.rstrip(file_ext)

    if file_ext == ".gz":
        with gzip.open(file_path, 'rb') as f_in:
            with open(unzip_path, 'wb') as f_out:
                f_out.write(f_in.read())
    elif file_ext == ".tar.gz":
        with tarfile.open(file_path, 'r:gz') as tar:
            tar.extractall(path=unzip_path)
    elif file_ext == ".zip":
        with zipfile.ZipFile(file_path, 'r') as zip_ref:
            zip_ref.extractall(unzip_path)
    elif file_ext == ".bz2":
        with bz2.open(file_path, 'rb') as f_in:
            with open(unzip_path, 'wb') as f_out:
                f_out.write(f_in.read())
    else:
        # Check compressed file and extract it
        unzip_path = os.path.dirname(os.path.abspath(file_path))

    return unzip_path

def create_index(input_seq_file):
    index = {}
    with open(input_seq_file, 'r') as f_seq:
        header = None
        sequence = ""
        for line in f_seq:
            if line.startswith(">"):
                if header:
                    index[header] = sequence
                header = line.strip().lstrip(">")
                sequence = ""
            else:
                sequence += line.strip()
        if header:
            index[header] = sequence
    return index

def find_duplicates(index):
    header_count = {}
    for header, _ in index.items():
        header_count[header] = header_count.get(header, 0) + 1
    return header_count

def remove_duplicates(index, header_count):
    unique_headers = set()
    removed_headers = set()

    for header, _ in index.items():
        if header_count[header] == 1:
            unique_headers.add(header)
        else:
            removed_headers.add(header)

    return unique_headers, removed_headers

def write_output_files(output_folder, index, unique_headers):
    outfasta_path = os.path.join(output_folder, "out_rmv_dupl.fasta")
    with open(outfasta_path, 'w') as f_outfasta:
        for header, sequence in index.items():
            if header in unique_headers:
                f_outfasta.write(f">{header}\n")
                f_outfasta.write(sequence + "\n")

def main():
    args = parse_args()
    logging.info(f"Input FASTA file: {args.input_seq}")
    logging.info(f"Output folder: {args.output}")
    logging.info(f"Number of CPUs: {args.t}")
    logging.info(f"Memory (GB): {args.mem}")

    # Read multi-line FASTA and convert to single-line
    sequences = read_fasta(args.input_seq)

    # Write the header counts to a file
    write_header_counts(sequences, args.output)

    # Remove duplicates
    unzip_path = unzip_file(args.input_seq)
    input_seq_file = os.path.join(unzip_path, os.path.basename(args.input_seq))

    index = create_index(input_seq_file)
    header_count = find_duplicates(index)
    unique_headers, _ = remove_duplicates(index, header_count)

    # Write the output fasta file
    if args.output:
        if not os.path.exists(args.output):
            os.makedirs(args.output)  # Create the output directory if it doesn't exist
    else:
        args.output = os.path.join(os.getcwd(), "output_files")

    write_output_files(args.output, index, unique_headers)

if __name__ == "__main__":
    main()

