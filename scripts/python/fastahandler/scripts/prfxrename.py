# FastaHandler rename selected and unselected IDs created by Hyungtaek Jung
# A multi-fasta (multiline) file to a single-line fasta with changing ID/Prefix name files (No space after the last IDs)
# Example command: python3 prexrename.py --input-seq test.fa --out test_out.fa --input-id new_id.txt (w/ optional for --t cpu and --mem memory)

#!/usr/bin/env python3
import argparse
import os
import multiprocessing as mp
import gzip
import bz2
import zipfile
import re
from pathlib import Path

# Function to decompress files based on extension
def decompress_file(file_path):
    if file_path.endswith('.gz'):
        opener = gzip.open
    elif file_path.endswith('.bz2'):
        opener = bz2.open
    elif file_path.endswith('.zip'):
        with zipfile.ZipFile(file_path, 'r') as zip_ref:
            zip_ref.extractall(Path(file_path).parent)
        return None
    else:
        opener = open

    with opener(file_path, 'rt') as f:
        return f.read()

# Function to process FASTA files
def process_fasta(content, id_map):
    fasta_header_pattern = re.compile(r'^>(.*?)\s*$')
    sequences = []
    current_sequence = []
    current_header = None

    for line in content.split('\n'):
        header_match = fasta_header_pattern.match(line)
        if header_match:
            if current_header:
                sequences.append((current_header, ''.join(current_sequence)))
            current_header = id_map.get(header_match.group(1), header_match.group(1))
            current_sequence = []
        else:
            current_sequence.append(line.upper())  # Case insensitive for sequence

    if current_header:
        sequences.append((current_header, ''.join(current_sequence)))

    return sequences

# Function to write output
def write_output(sequences, output_file):
    with open(output_file, 'w') as f:
        for header, seq in sequences:
            f.write(f">{header}\n{seq}\n")

# Function to read ID map
def read_id_map(file_path):
    id_map = {}
    with open(file_path) as f:
        for line in f:
            original_id, new_id = line.strip().split('\t')
            id_map[original_id] = new_id
    return id_map

# Main function
def main():
    parser = argparse.ArgumentParser(description='prfxrename - FASTA file processing and renaming tool')
    parser.add_argument('--input-seq', required=True, help='Path to the input sequence file')
    parser.add_argument('--input-id', required=True, help='Path to the file with new IDs')
    parser.add_argument('--out', required=True, help='Path to the output file')
    parser.add_argument('--t', type=int, default=1, help='Number of CPUs (default: 1)')
    parser.add_argument('--mem', type=int, default=10, help='Memory in GB (default: 10)')
    args = parser.parse_args()

    # Validate file extensions
    if not (args.input_seq.endswith(('.fasta', '.fa', '.fas', '.fna', '.tar.gz', '.gz', '.zip', '.bz2'))):
        raise ValueError("Not a proper file format. Please provide a readable file and sequence format to call the prfxrename module.")

    # Decompress file if needed
    if args.input_seq.endswith(('.tar.gz', '.gz', '.zip', '.bz2')):
        content = decompress_file(args.input_seq)
        if content is None:  # Special case for zip files
            args.input_seq = args.input_seq.replace('.zip', '')
            with open(args.input_seq, 'r') as f:
                content = f.read()
    else:
        with open(args.input_seq, 'r') as f:
            content = f.read()

    # Read ID map
    id_map = read_id_map(args.input_id)

    # Process FASTA
    sequences = process_fasta(content, id_map)

    # Write output
    write_output(sequences, args.out)

    print("Processing completed.")

if __name__ == "__main__":
    main()

