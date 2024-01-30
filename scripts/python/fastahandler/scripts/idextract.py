# FastaHandler extract IDs created by Hyungtaek Jung
# Extract matched ID/header sequences from a multi-fasta (multiline) file to be written in a single-line fasta
# Must create the input header file after removing ">"
# Example command: python3 idextract.py --input-seq test.fa --input-hdr test_id.txt --out test_out.fa (w/ optional for --t cpu and --mem memory)


#!/usr/bin/env python3
import argparse
import os
import sys
import gzip
import tarfile
import zipfile
import bz2
import logging

logging.basicConfig(filename="idextract.log", level=logging.INFO, format="%(asctime)s - %(message)s")

def parse_args():
    parser = argparse.ArgumentParser(description="Extract sequences based on headers from a multi-line FASTA file.")
    parser.add_argument("--input-seq", required=True, help="Indicate the input FASTA file and path.")
    parser.add_argument("--input-hdr", required=True, help="Indicate the input header file and path.")
    parser.add_argument("--out", required=False, help="Indicate the output FASTA file and path.")
    parser.add_argument("--t", type=int, default=1, help="Number of CPUs with only numbers (integers). Default: 1.")
    parser.add_argument("--mem", type=int, default=10, help="Number of memory in gigabytes with only numbers (integers). Default: 10.")
    args = parser.parse_args()

    if not os.path.exists(args.input_seq):
        sys.exit("Error: Input FASTA file not found or inaccessible.")

    if not os.path.exists(args.input_hdr):
        sys.exit("Error: Input header file not found or inaccessible or containing '>'.")

    if args.out:
        out_dir = os.path.dirname(os.path.abspath(args.out))
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)  # Create the output directory if it doesn't exist
    else:
        args.out = os.path.join(os.getcwd(), "output_files.fasta")

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

def extract_sequences(input_seq_file, header_set, output_file):
    with open(input_seq_file, 'r') as f_seq, open(output_file, 'w') as f_out:
        header = None
        sequence = ""
        for line in f_seq:
            if line.startswith(">"):
                if header and header in header_set:
                    f_out.write(f">{header}\n")
                    f_out.write(sequence.replace("\n", "") + "\n")
                header = line.strip().lstrip(">")
                sequence = ""
            else:
                sequence += line.strip()
        if header and header in header_set:
            f_out.write(f">{header}\n")
            f_out.write(sequence.replace("\n", "") + "\n")

def main():
    args = parse_args()
    logging.info(f"Input FASTA file: {args.input_seq}")
    logging.info(f"Input header file: {args.input_hdr}")
    logging.info(f"Output FASTA file: {args.out}")
    logging.info(f"Number of CPUs: {args.t}")
    logging.info(f"Memory (GB): {args.mem}")

    unzip_path_seq = unzip_file(args.input_seq)
    input_seq_file = os.path.join(unzip_path_seq, os.path.basename(args.input_seq))

    unzip_path_header = unzip_file(args.input_hdr)
    input_header_file = os.path.join(unzip_path_header, os.path.basename(args.input_hdr))

    with open(input_header_file, 'r') as f_header:
        header_set = set(line.strip() for line in f_header if line.strip())

    extract_sequences(input_seq_file, header_set, args.out)

if __name__ == "__main__":
    main()

