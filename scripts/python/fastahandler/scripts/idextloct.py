# FastaHandler extract a single ID and location created by Hyungtaek Jung
# Extract matched ID/header sequences and their locations from a multi-fasta (multiline) file to be written in a single-line fasta
# Must create the input header file after removing ">"
# Note, users must provide -1 less numbers for start and end location. 
# Example command: python3 idextloct.py --input-seq test.fa --hdr-id Seq%_3 --start 30 --end 40 --out test_out.fa (w/ optional for --t cpu and --mem memory)


#!/usr/bin/env python3
import argparse
import os
import sys
import gzip
import tarfile
import zipfile
import bz2
import logging

logging.basicConfig(filename="idextloc.log", level=logging.INFO, format="%(asctime)s - %(message)s")

def parse_args():
    parser = argparse.ArgumentParser(description="Extract sequences based on headers and location from a multi-line FASTA file.")
    parser.add_argument("--input-seq", required=True, help="Indicate the input FASTA file and path.")
    parser.add_argument("--hdr-id", required=True, help="Indicate the header ID to extract without '>'")
    parser.add_argument("--start", type=int, required=True, help="Indicate the start location to extract (1-based).")
    parser.add_argument("--end", type=int, required=True, help="Indicate the end location to extract (1-based).")
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
        # Check the comressed file and extract it
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

def main():
    args = parse_args()
    logging.info(f"Input FASTA file: {args.input_seq}")
    logging.info(f"Header ID: {args.hdr_id}")
    logging.info(f"Start location: {args.start}")
    logging.info(f"End location: {args.end}")
    logging.info(f"Output FASTA file: {args.out}")
    logging.info(f"Number of CPUs: {args.t}")
    logging.info(f"Memory (GB): {args.mem}")

    unzip_path = unzip_file(args.input_seq)
    input_seq_file = os.path.join(unzip_path, os.path.basename(args.input_seq))

    index = create_index(input_seq_file)

    if args.hdr_id not in index:
        sys.exit("Error: Header ID does not match the indexed sequence.")

    sequence = index[args.hdr_id]
    sequence_length = len(sequence)
    start = args.start - 1
    end = args.end - 1

    if start >= sequence_length or end >= sequence_length or start > end:
        sys.exit("Error: Start and/or end locations do not match the indexed sequence.")

    extracted_sequence = sequence[start:end+1]

    if args.out:
        out_dir = os.path.dirname(os.path.abspath(args.out))
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)  # Create the output directory if it doesn't exist
    else:
        args.out = os.path.join(os.getcwd(), "output_files.fasta")

    with open(args.out, 'w') as f_out:
        f_out.write(f">{args.hdr_id}\n")
        f_out.write(extracted_sequence + "\n")

if __name__ == "__main__":
    main()

