# FastaHandler reverse complement created by Hyungtaek Jung
# Make reverse complement sequences from a multi-fasta (multiline) file to be written in a single-line fasta
# Example command: python3 revcomplt.py --input-seq test.fa --out test_out.fa (w/ optional for --t cpu and --mem memory)


#!/usr/bin/env python3
import argparse
import os
import sys
import gzip
import tarfile
import zipfile
import bz2
import logging

logging.basicConfig(filename="revcomplt.log", level=logging.INFO, format="%(asctime)s - %(message)s")

def parse_args():
    parser = argparse.ArgumentParser(description="Generate the reverse complement of sequences from a multi-line FASTA file.")
    parser.add_argument("--input-seq", required=True, help="Indicate the input FASTA file and path.")
    parser.add_argument("--out", required=False, help="Indicate the output FASTA file and path.")
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
        # If the file is not in a supported compressed format, assume it's already uncompressed
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

def reverse_complement(sequence):
    complement = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N",
                  "a": "t", "t": "a", "c": "g", "g": "c", "n": "n"}
    rev_complement = "".join(complement[base] for base in reversed(sequence))
    return rev_complement


def main():
    args = parse_args()
    logging.info(f"Input FASTA file: {args.input_seq}")
    logging.info(f"Output FASTA file: {args.out}")
    logging.info(f"Number of CPUs: {args.t}")
    logging.info(f"Memory (GB): {args.mem}")

    unzip_path = unzip_file(args.input_seq)
    input_seq_file = os.path.join(unzip_path, os.path.basename(args.input_seq))

    index = create_index(input_seq_file)

    if args.out:
        out_dir = os.path.dirname(os.path.abspath(args.out))
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)  # Create the output directory if it doesn't exist
    else:
        args.out = os.path.join(os.getcwd(), "output_revcomplement.fasta")

    with open(args.out, 'w') as f_out:
        for header, sequence in index.items():
            rev_complement_sequence = reverse_complement(sequence)
            f_out.write(f">{header}\n")
            f_out.write(rev_complement_sequence + "\n")

if __name__ == "__main__":
    main()

