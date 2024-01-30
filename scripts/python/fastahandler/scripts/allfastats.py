# FastaHandler a multi-line fasta created by Hyungtaek Jung
# A multi-fasta (multiline) file to calculate assembly stats
# Example command: python3 allfastats.py --input-seq test1.fa --out all_stats_out.txt (w/ optional for --t cpu and --mem memory)


#!/usr/bin/env python3

import argparse
import os
import re
import gzip
import logging
from collections import Counter
import numpy as np
import pandas as pd

# Set up logging
logging.basicConfig(filename='allfastats.log', level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')

def unzip_file(input_file, output_dir):
    # Check compressed file and extract it
    if input_file.lower().endswith(('.gz', '.tar.gz', '.zip', '.bz2')):
        with gzip.open(input_file, 'rb') as f_in:
            with open(os.path.join(output_dir, os.path.basename(input_file)[:-3]), 'wb') as f_out:
                f_out.write(f_in.read())
        return os.path.join(output_dir, os.path.basename(input_file)[:-3])
    return input_file

def parse_fasta(file_path):
    # Check compressed file and extract it
    file_path = unzip_file(file_path, os.path.dirname(file_path))

    sequences = {}
    current_seq_id = None
    with open(file_path, "r") as file:
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                current_seq_id = line[1:]
                sequences[current_seq_id] = ""
            else:
                sequences[current_seq_id] += line

    # Convert multiline sequences to single-line format
    for seq_id, sequence in sequences.items():
        sequences[seq_id] = sequence.replace(" ", "").upper()

    return sequences

def calculate_summary_stats(file_path):
    sequences = parse_fasta(file_path)
    genome = "".join(sequences.values())
    
    input_id = os.path.basename(file_path)
    assembled_genome = len(genome)
    contigs = len(sequences)
    gap_count = genome.count("-")
    uncertain_bp = sum(1 for char in genome if char not in "ATGCN-")

    # Count occurrences of each nucleotide
    nucleotide_counts = Counter(genome)
    total_bp = len(genome)

    # Calculate percentages
    gc_percent = (nucleotide_counts["G"] + nucleotide_counts["C"]) / total_bp * 100
    a_percent = nucleotide_counts["A"] / total_bp * 100
    t_percent = nucleotide_counts["T"] / total_bp * 100
    g_percent = nucleotide_counts["G"] / total_bp * 100
    c_percent = nucleotide_counts["C"] / total_bp * 100
    n_percent = nucleotide_counts["N"] / total_bp * 100

    # Calculate N50
    sorted_contigs = sorted(sequences.values(), key=len, reverse=True)
    total_length = sum(len(contig) for contig in sorted_contigs)
    half_total_length = total_length / 2
    n50 = None
    for contig_length in map(len, sorted_contigs):
        half_total_length -= contig_length
        if half_total_length <= 0:
            n50 = contig_length
            break

    # Calculate sequence length ranges
    range_counts = Counter()
    for length in map(len, sequences.values()):
        if length < 200:
            range_counts["200bp"] += 1
        elif 200 < length < 1000:
            range_counts["1Kbp"] += 1
        elif 1000 < length < 5000:
            range_counts["5Kbp"] += 1
        elif 5000 < length < 10000:
            range_counts["10Kbp"] += 1
        elif 10000 < length < 50000:
            range_counts["50Kbp"] += 1
        elif 50000 < length < 100000:
            range_counts["100Kbp"] += 1
        elif 100000 < length < 500000:
            range_counts["500Kbp"] += 1
        elif 500000 < length < 1000000:
            range_counts["1Mbp"] += 1
        elif 1000000 < length < 5000000:
            range_counts["5Mbp"] += 1
        else:
            range_counts["Over10Mbp"] += 1

    logging.info("Summary statistics calculated successfully.")
    return {
        "InputID": input_id,
        "AssembledGenome": assembled_genome,
        "Contigs": contigs,
        "GC(%)": round(gc_percent, 2),
        "A(%)": round(a_percent, 2),
        "T(%)": round(t_percent, 2),
        "G(%)": round(g_percent, 2),
        "C(%)": round(c_percent, 2),
        "N(%)": round(n_percent, 2),
        "Ncount": nucleotide_counts["N"],
        "N50": n50,
        "Max": len(sorted_contigs[0]),
        "Min": len(sorted_contigs[-1]),
        "GapCount(-)": gap_count,
        "Uncertain(bp)": uncertain_bp,
        "200bp": range_counts["200bp"],
        "1Kbp": range_counts["1Kbp"],
        "5Kbp": range_counts["5Kbp"],
        "10Kbp": range_counts["10Kbp"],
        "50Kbp": range_counts["50Kbp"],
        "100Kbp": range_counts["100Kbp"],
        "500Kbp": range_counts["500Kbp"],
        "1Mbp": range_counts["1Mbp"],
        "5Mbp": range_counts["5Mbp"],
        "Over10Mbp": range_counts["Over10Mbp"]
    }

def generate_summary_output(summary_stats, output_file):
    headers = [
        "InputID", "AssembledGenome", "Contigs", "GC(%)", "A(%)", "T(%)", "G(%)", "C(%)", "N(%)", "Ncount", "N50",
        "Max", "Min", "GapCount(-)", "Uncertain(bp)", "200bp", "1Kbp", "5Kbp", "10Kbp", "50Kbp", "100Kbp", "500Kbp",
        "1Mbp", "5Mbp", "Over10Mbp"
    ]

    data = [[summary_stats[header] for header in headers]]
    df = pd.DataFrame(data, columns=headers)
    df.to_csv(output_file, sep='\t', index=False)

    logging.info(f"Summary output written to {output_file} successfully.")

def main():
    parser = argparse.ArgumentParser(description="Comprehensive Python script for analyzing fasta files.")
    parser.add_argument("--input-seq", required=True, help="Input fasta file path.")
    parser.add_argument("--out", required=True, help="Output file path for summary results.")
    parser.add_argument("--t", type=int, default=1, help="Number of CPUs.")
    parser.add_argument("--mem", type=int, default=10, help="Amount of memory in gigabytes.")

    args = parser.parse_args()

    # Check for valid CPU and memory values
    if args.t < 1 or args.mem < 1:
        print("Not an integer number. Please provide only readable integer numbers to call the proper CPUs and Memory.")
        return
    logging.info(f"Using {args.t} CPUs and {args.mem}GB of memory.")

    # Check and unzip input file if necessary
    if not os.path.isfile(args.input_seq):
        raise FileNotFoundError("Input file not found.")
    logging.info(f"Input file: {args.input_seq}")

    # Calculate summary statistics
    summary_stats = calculate_summary_stats(args.input_seq)

    # Generate log file
    logging.info(f"Summary statistics: {summary_stats}")

    # Generate summary output
    generate_summary_output(summary_stats, args.out)

if __name__ == "__main__":
    main()
