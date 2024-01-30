# Open Reading Frames / Transcription May Begin Anywhere
# 6 possible reading frames includng forward & reverse complement sequences
# ">>>" and "<<<" tags to highlight the longest_nt_seq and longest_protein_seq
# Fasta file header must be no "space" and convert an automatic singleline
# Need to install Biopython
# A total of four outputs (Edited from 6f4.py)
# Example usage: python script_CC.py --input-seq test_files.fasta.gz --num-processes 4


#!/usr/bin/env python3

import os
import argparse
import glob
import gzip
import bz2
import tarfile
import zipfile
from Bio import SeqIO
from Bio.Seq import Seq
from multiprocessing import Pool
from functools import partial

# Function to extract and read sequences from compressed files
def extract_and_read_sequences(file_path):
    if file_path.endswith(".gz"):
        opener = gzip.open
    elif file_path.endswith(".bz2"):
        opener = bz2.open
    elif file_path.endswith(".zip"):
        with zipfile.ZipFile(file_path, "r") as z:
            z.extractall(os.path.dirname(file_path))
            file_list = z.namelist()
            return [SeqIO.parse(os.path.join(os.path.dirname(file_path), fname), "fasta") for fname in file_list]
    elif file_path.endswith(".tar.gz"):
        with tarfile.open(file_path, "r:gz") as tar:
            tar.extractall(os.path.dirname(file_path))
            return [SeqIO.parse(os.path.join(os.path.dirname(file_path), member.name), "fasta") for member in tar.getmembers() if member.isfile()]
    elif file_path.endswith(('.fasta', '.fa', '.fas', '.fna')):
        return [SeqIO.parse(file_path, "fasta")]
    else:
        raise ValueError("Not a proper sequence file format. Please provide a readable file and sequence format.")

# Function to find ORFs in a sequence
def find_orfs(sequence):
    start_codon = "ATG"
    stop_codons = ["TAA", "TAG", "TGA"]
    orfs = []

    forward_sequence = sequence
    reverse_sequence = forward_sequence.reverse_complement()

    for seq in [forward_sequence, reverse_sequence]:
        for i in range(0, len(seq), 3):
            if seq[i:i + 3] == start_codon:
                for j in range(i + 3, len(seq), 3):
                    if seq[j:j + 3] in stop_codons:
                        orfs.append(str(seq[i:j + 3]))
                        break

    return list(set(orfs))

# Function to process a single record
def process_record(record, output_dir):
    header = record.id
    sequence = record.seq

    try:
        orfs_nt = find_orfs(sequence)
        orfs_protein = [str(Seq(orf).translate(to_stop=True)) for orf in orfs_nt]

        longest_nt_seq = max(orfs_nt, key=len) if orfs_nt else None
        longest_protein_seq = max(orfs_protein, key=len) if orfs_protein else None

        output_files = {
            "nt": os.path.join(output_dir, "out_nt.fasta"),
            "prot": os.path.join(output_dir, "out_prot.fasta"),
            "nt_sltd": os.path.join(output_dir, "out_nt_sltd.fasta"),
            "prot_sltd": os.path.join(output_dir, "out_prot_sltd.fasta")
        }

        for orf_type, orfs in [("nt", orfs_nt), ("prot", orfs_protein)]:
            with open(output_files[orf_type], "a") as f:
                f.write(f">{header}\n" + "\n".join(orfs) + "\n")

        if longest_nt_seq:
            with open(output_files["nt_sltd"], "a") as f:
                f.write(f">{header}_longest_nt\n{longest_nt_seq}\n")

        if longest_protein_seq:
            with open(output_files["prot_sltd"], "a") as f:
                f.write(f">{header}_longest_protein\n{longest_protein_seq}\n")

    except Exception as e:
        print(f"Error processing record with header: {header}. Error: {e}")

# Function to process files
def process_files(input_files, output_dir, threads):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for input_file in input_files:
        try:
            records_list = extract_and_read_sequences(input_file)
            for records in records_list:
                with Pool(processes=threads) as pool:
                    process_func = partial(process_record, output_dir=output_dir)
                    pool.map(process_func, records)
        except Exception as e:
            raise ValueError(f"Error reading {input_file}. Ensure it's a proper FASTA file. Error: {e}")

# Main function to handle input and output
def main(input_seq, output_dir, threads):
    input_files = glob.glob(input_seq)
    if not input_files:
        raise ValueError("No files found with the given pattern.")

    process_files(input_files, output_dir, threads)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process multi-line FASTA files.")
    parser.add_argument("--input-seq", required=True, help="Path or pattern to match input FASTA files (e.g., '/path/*.fasta').")
    parser.add_argument("--out", required=True, help="Output directory for all outcome files.")
    parser.add_argument("--t", type=int, default=1, help="Number of CPU threads (default: 1).")
    parser.add_argument("--mem", type=int, default=10, help="Memory in GB (default: 10).")
    args = parser.parse_args()

    if args.t < 1:
        raise ValueError("Please provide a positive integer number for CPU threads.")

    # Note: The memory argument is accepted but not actively used in memory management in this script.
    main(args.input_seq, args.out, args.t)

