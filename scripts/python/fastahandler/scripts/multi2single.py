# FastaHandler multiline to single-line created by Hyungtaek Jung
# A multi-fasta (multiline) file to a single-line fasta
# Example command: python3 multi2single.py --input-seq test.fa --out test_out.fa (w/ optional for --t cpu and --mem memory)


#!/usr/bin/env python3
import argparse
import os
import sys
import gzip
import tarfile
import zipfile
import bz2
import logging

logging.basicConfig(filename="multi2single.log", level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

def parse_args():
    parser = argparse.ArgumentParser(description="Convert multi-line FASTA to single-line FASTA.")
    parser.add_argument("--input-seq", required=True, help="Input file path.")
    parser.add_argument("--out", required=False, help="Output file path.")
    parser.add_argument("--t", type=int, default=1, help="Number of CPUs. Default: 1.")
    parser.add_argument("--mem", type=int, default=10, help="Memory in gigabytes. Default: 10.")
    args = parser.parse_args()

    if not os.path.exists(args.input_seq):
        sys.exit("Error: Input file not found or inaccessible.")

    if args.out:
        out_dir = os.path.dirname(os.path.abspath(args.out))
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)  # Create the output directory if it doesn't exist
    else:
        args.out = os.path.join(os.getcwd(), "renameid_seq_T1.fa")

    if not isinstance(args.t, int) or args.t <= 0:
        sys.exit("Error: Number of CPUs must be a positive integer.")

    if not isinstance(args.mem, int) or args.mem <= 0:
        sys.exit("Error: Memory must be a positive integer.")

    return args

def unzip_file(file_path):
    file_ext = os.path.splitext(file_path)[-1]
    unzip_path = file_path.rstrip(file_ext)

    if file_ext == ".gz":
        with gzip.open(file_path, 'rt') as f_in:
            with open(unzip_path, 'w') as f_out:
                f_out.write(f_in.read())
    elif file_ext == ".tar.gz":
        with tarfile.open(file_path, 'r:gz') as tar:
            tar.extractall(path=unzip_path)
    elif file_ext == ".zip":
        with zipfile.ZipFile(file_path, 'r') as zip_ref:
            zip_ref.extractall(unzip_path)
    elif file_ext == ".bz2":
        with bz2.open(file_path, 'rt') as f_in:
            with open(unzip_path, 'w') as f_out:
                f_out.write(f_in.read())
    else:
        unzip_path = os.path.dirname(os.path.abspath(file_path))

    return unzip_path

def multi2singleline(input_file, output_file):
    with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
        header = None
        sequence = ""
        for line in f_in:
            if line.startswith(">"):
                if header is not None:
                    f_out.write(header + "\n")
                    f_out.write(sequence + "\n")
                header = line.strip()
                sequence = ""
            else:
                sequence += line.strip()
        if header is not None and sequence:
            f_out.write(header + "\n")
            f_out.write(sequence + "\n")

def main():
    args = parse_args()
    logging.info(f"Input file: {args.input_seq}")
    logging.info(f"Output file: {args.out}")
    logging.info(f"Number of CPUs: {args.t}")
    logging.info(f"Memory (GB): {args.mem}")

    # Convert to a single line
    unzip_path = unzip_file(args.input_seq)
    single_line_output = args.out or os.path.join(os.getcwd(), "output_singleline.fasta")
    multi2singleline(args.input_seq, single_line_output)

if __name__ == "__main__":
    main()
