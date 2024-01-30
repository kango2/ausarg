# FASTAhandler extract multiple IDs and locations created by Hyungtaek Jung
# Extract multiple matched ID/header sequences and their locations from a multi-fasta (multiline) file to be written in a single-line fasta
# Must create the input header file after removing ">"
# Must provide with a tap-separated ID/header and their locations (start/end)
# Note, users must provide -1 less numbers for start and end location. 
# Example command: python3 idextloctmlt.py --input-seq test.fa --input-ext hdrmulti_id.txt --out test_out.fa (w/ optional for --t cpu and --mem memory)


#!/usr/bin/env python3
import argparse
import os
import gzip
import tarfile
import zipfile
import bz2
import logging
import re
import sys

logging.basicConfig(filename="idextractlocamulti.log", level=logging.INFO, format="%(asctime)s - %(message)s")

def parse_args():
    parser = argparse.ArgumentParser(description="Extract sequences based on headers and location from a multi-line FASTA file using multiple header and location information.")
    parser.add_argument("--input-seq", required=True, help="Indicate the input FASTA file and path.")
    parser.add_argument("--input-ext", required=True, help="Indicate the input file (text, tab-separated) and path for header, start, and end locations.")
    parser.add_argument("--out", required=False, help="Indicate the output FASTA file and path.")
    parser.add_argument("--t", type=int, default=1, help="Number of CPUs with only numbers (integers). Default: 1.")
    parser.add_argument("--mem", type=int, default=10, help="Number of memory in gigabytes with only numbers (integers). Default: 10.")
    args = parser.parse_args()

    if not os.path.exists(args.input_seq):
        sys.exit("Error: Input FASTA file not found or inaccessible.")

    if not os.path.exists(args.input_ext):
        sys.exit("Error: Input extract file not found or inaccessible.")

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
    current_header = None
    current_sequence = []

    with open(input_seq_file, 'r') as f_seq:
        for line in f_seq:
            line = line.strip()
            if line.startswith(">"):
                if current_header is not None:
                    index[current_header] = "".join(current_sequence)
                current_header = re.sub(r'^>', '', line)  # Remove '>' and any leading/trailing whitespaces
                current_sequence = []
            else:
                current_sequence.append(line)

        if current_header is not None:
            index[current_header] = "".join(current_sequence)

    return index

def read_extract_file(input_extract_file):
    extract_info = []
    with open(input_extract_file, 'r') as f_extract:
        for line_num, line in enumerate(f_extract, 1):
            fields = line.strip().split("\t")
            if len(fields) != 3:
                logging.warning(f"Invalid line {line_num} in the extract file. Skipping.")
                continue

            header, start, end = fields

            try:
                start = int(start)
                end = int(end)
            except ValueError:
                logging.warning(f"Invalid start/end position in line {line_num} of the extract file. Skipping extraction for header '{header}'.")
                continue

            extract_info.append((header, start, end))
    return extract_info

def main():
    args = parse_args()
    logging.info(f"Input FASTA file: {args.input_seq}")
    logging.info(f"Input extract file: {args.input_ext}")
    logging.info(f"Output FASTA file: {args.out}")
    logging.info(f"Number of CPUs: {args.t}")
    logging.info(f"Memory (GB): {args.mem}")

    # Handle decompression of input FASTA file (if needed)
    unzip_path = unzip_file(args.input_seq)
    input_seq_file = os.path.join(unzip_path, os.path.basename(args.input_seq))

    # Create an index for the input FASTA file
    index = create_index(input_seq_file)

    # Read the extract file containing header, start, and end locations
    extract_info = read_extract_file(args.input_ext)

    # Perform sequence extraction and write to the output file
    with open(args.out, 'w') as f_out:
        for header, start, end in extract_info:
            # Check if the header is present in the index
            if header not in index:
                logging.warning(f"Header ID '{header}' not found in the indexed sequence.")
                continue

            # Get the full sequence for the header
            sequence = index[header]
            sequence_length = len(sequence)

            # Check if start and end locations are valid
            if start >= sequence_length or end >= sequence_length or start > end:
                logging.warning(f"Invalid start and/or end locations for header '{header}'. Skipping extraction.")
                continue

            # Extract the sequence based on the start and end locations
            extracted_sequence = sequence[start:end+1]

            # Write the extracted sequence to the output file in FASTA format
            f_out.write(f">{header}\n")
            f_out.write(extracted_sequence + "\n")

if __name__ == "__main__":
    main()

