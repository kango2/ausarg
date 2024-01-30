# FASTAhandler find only duplicated ID created by Hyungtaek Jung
# Find and count the duplicate ID and sequence from a multi-fasta (multiline) file
# Example command: python3 findcntdupl.py --input-seq test.fa --out test_out.txt (w/ optional for --t cpu and --mem memory)


#!/usr/bin/env python3

import argparse
import os
import gzip
import tarfile
import zipfile
import bz2
import logging

logging.basicConfig(filename="findcntdupl.log", level=logging.INFO, format="%(asctime)s - %(message)s")

def parse_args():
    parser = argparse.ArgumentParser(description="Find and count duplicated headers in a multi-line FASTA file.")
    parser.add_argument("--input-seq", required=True, help="Indicate the input FASTA file and path.")
    parser.add_argument("--out", required=False, help="Indicate the output text file and path for summary.")
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

    if file_ext in (".gz", ".tar.gz", ".bz2"):
        with get_compression_open(file_path) as f_in:
            with open(unzip_path, 'wb') as f_out:
                f_out.write(f_in.read())
    elif file_ext == ".zip":
        with zipfile.ZipFile(file_path, 'r') as zip_ref:
            zip_ref.extractall(unzip_path)
    else:
        # Check compressed file and extract it
        unzip_path = os.path.dirname(os.path.abspath(file_path))

    return unzip_path

def get_compression_open(file_path):
    if file_path.endswith(".gz"):
        return gzip.open(file_path, 'rb')
    elif file_path.endswith(".bz2"):
        return bz2.open(file_path, 'rb')
    elif file_path.endswith(".tar.gz"):
        return tarfile.open(file_path, 'r:gz')
    else:
        raise ValueError(f"Unsupported compression format for file: {file_path}")

def create_index(input_seq_file):
    index = {}
    with open(input_seq_file, 'r') as f_seq:
        header = None
        for line in f_seq:
            if line.startswith(">"):
                header = line.strip()[1:].strip()  # Remove ">" and leading/trailing whitespaces
                index[header] = index.get(header, 0) + 1
            # Counting each sequence line is not necessary, as headers are unique
    return index

def main():
    args = parse_args()
    logging.info(f"Input FASTA file: {args.input_seq}")
    logging.info(f"Output summary file: {args.out}")
    logging.info(f"Number of CPUs: {args.t}")
    logging.info(f"Memory in GB: {args.mem}")

    unzip_path = unzip_file(args.input_seq)
    input_seq_file = os.path.join(unzip_path, os.path.basename(args.input_seq))

    index = create_index(input_seq_file)

    if args.out:
        out_dir = os.path.dirname(os.path.abspath(args.out))
        os.makedirs(out_dir, exist_ok=True)  # Create the output directory if it doesn't exist
    else:
        args.out = os.path.join(os.getcwd(), "output_files.txt")

    with open(args.out, 'w') as f_out:
        f_out.write("Headers\tCounts\n")
        for header, count in index.items():
            if count > 1:
                f_out.write(f"{header}\t{count}\n")

if __name__ == "__main__":
    main()

