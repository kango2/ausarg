# FastaHandler multiline fasta assembly stats PyFastaStats created by Hyungtaek Jung
# Unlimited multiple-multi-fasta (multiline) files to calculate assembly stats
# Example command: python3 multiasmstats_unlimt.py --input-seqs test1.fa test2.fa test3.fa test4.fa test5.fa --out multi_stats_ulmt.txt (w/ optional for --t cpu and --mem memory)


#!/usr/bin/env python3

import argparse
import os
import re
import gzip
import logging
from collections import Counter
import pandas as pd
from concurrent.futures import ProcessPoolExecutor

# Set up logging
logging.basicConfig(filename='multiasmstats.log', level=logging.INFO,
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

    # Convert multiline sequences to single-line format and uppercase
    sequences = {seq_id: sequence.replace(" ", "").upper() for seq_id, sequence in sequences.items()}
    return sequences

def calculate_summary_stats(file_path):
    sequences = parse_fasta(file_path)
    genome = "".join(sequences.values())

    # Calculate basic statistics
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
    ranges = [200, 1000, 5000, 10000, 50000, 100000, 500000, 1000000, 5000000]
    range_counts = {f"{r}bp": 0 for r in ranges}
    range_counts["Over10Mbp"] = 0

    for length in map(len, sequences.values()):
        for r in ranges:
            if length < r:
                range_counts[f"{r}bp"] += 1
                break
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
        **range_counts  # include range counts dynamically
    }

def generate_summary_output(summary_stats_list, output_file):
    # Convert the summary stats to a DataFrame
    df = pd.DataFrame(summary_stats_list).set_index("InputID")

    # Transpose the DataFrame to get the desired format
    df = df.transpose()

    # Write DataFrame to a tab-separated file
    df.to_csv(output_file, sep='\t')

    logging.info(f"Summary output written to {output_file} successfully.")

def main():
    parser = argparse.ArgumentParser(description="Comprehensive Python script for analyzing fasta files.")
    parser.add_argument("--input-seqs", nargs='+', required=True, help="Input fasta file paths for sequences.")
    parser.add_argument("--out", required=True, help="Output file path for summary results.")
    parser.add_argument("--t", type=int, default=1, help="Number of CPUs.")
    parser.add_argument("--mem", type=int, default=10, help="Amount of memory in gigabytes.")

    args = parser.parse_args()

    # Check for valid CPU and memory values
    if args.t < 1 or args.mem < 1:
        print("Not an integer number. Please provide only readable integer numbers to call the proper CPUs and Memory.")
        return
    logging.info(f"Using {args.t} CPUs and {args.mem}GB of memory.")

    # Check and unzip input files if necessary
    input_files = args.input_seqs
    for file_path in input_files:
        if not os.path.isfile(file_path):
            raise FileNotFoundError(f"Input file not found: {file_path}")
        logging.info(f"Input file: {file_path}")

    # Calculate summary statistics for each input file using multiprocessing
    with ProcessPoolExecutor(max_workers=args.t) as executor:
        summary_stats_list = list(executor.map(calculate_summary_stats, input_files))

    # Generate log file
    logging.info("Summary statistics for each input file:")
    for file_path, summary_stats in zip(input_files, summary_stats_list):
        logging.info(f"{file_path}: {summary_stats}")

    # Generate summary output
    generate_summary_output(summary_stats_list, args.out)

if __name__ == "__main__":
    main()

