# FastaHandler subset and filter length created by Hyungtaek Jung
# Subset based on sequence length from a multi-fasta (multiline) file
# Example usage: python subsetfa.py --input-seq test_mRNA1.fasta --filter 50 --out output_subset.fasta (w/ optional for --t cpu and --mem memory)


import argparse
import logging
import os
import time
from multiprocessing import Pool, cpu_count
from Bio import SeqIO

def extract_sequences(file_path):
    if file_path.endswith(".gz"):
        with gzip.open(file_path, "rt") as f:
            return f.read()
    elif file_path.endswith(".zip"):
        with zipfile.ZipFile(file_path, "r") as z:
            file_list = z.namelist()
            if len(file_list) == 1:  # Assuming only one file in the zip archive
                with z.open(file_list[0], "r") as f:
                    return f.read().decode()
    elif file_path.endswith(".bz2"):
        with bz2.open(file_path, "rt") as f:
            return f.read()
    else:
        raise ValueError("Not a supported compressed file format.")

def find_sequence_length(sequence):
    cleaned_seq = "".join(sequence.upper().split())  # Remove spaces and tabs
    return len(cleaned_seq)

def filter_and_write(args):
    record, filter_length, output_file = args
    if len(record.seq) >= filter_length:
        seq_str = str(record.seq)
        with open(output_file, "a") as f_out:
            f_out.write(f">{record.id}\n{seq_str}\n")

def main(input_seq, filter_length, output_file, num_processes, memory):
    # Create or truncate the output file
    with open(output_file, "w"):
        pass

    # Configure logging
    logging.basicConfig(filename="subsetfa_log.txt", level=logging.INFO,
                        format="%(asctime)s - %(levelname)s: %(message)s")
    logging.info("Running subsetfasta.py with the following parameters:")
    logging.info(f"Input file: {input_seq}")
    logging.info(f"Filter length: {filter_length}")
    logging.info(f"Output file: {output_file}")
    logging.info(f"Number of CPUs: {num_processes}")
    logging.info(f"Memory (GB): {memory}")

    # Read the input FASTA file and get the list of records
    try:
        records = list(SeqIO.parse(input_seq, "fasta"))
    except Exception as e:
        logging.error("Not a proper sequence file format. Please provide a readable FASTA file.")
        raise ValueError("Not a proper sequence file format. Please provide a readable FASTA file.") from e

    # Process records using multiprocessing
    if num_processes is None or num_processes <= 0:
        num_processes = cpu_count()

    args_list = [(record, filter_length, output_file) for record in records]
    with Pool(processes=num_processes) as pool:
        pool.map(filter_and_write, args_list)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Subset FASTA files based on sequence length.")
    parser.add_argument("--input-seq", required=True, help="Path to the input sequence file.")
    parser.add_argument("--filter", type=int, required=True, help="Filter sequences with length greater than or equal to this value.")
    parser.add_argument("--out", required=True, help="Output file for the filtered sequences in FASTA format.")
    parser.add_argument("--t", type=int, default=1, help="Number of CPUs for performance (default: 1).")
    parser.add_argument("--mem", type=int, default=10, help="Memory in gigabytes (default: 10).")
    args = parser.parse_args()

    start_time = time.time()
    main(args.input_seq, args.filter, args.out, args.t, args.mem)
    end_time = time.time()

    print(f"Script executed in {end_time - start_time:.2f} seconds.")

