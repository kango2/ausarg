import os
import sys
import hashlib
import gzip
import csv
from Bio import SeqIO

def md5(file_path, is_compressed=False):
    hash_md5 = hashlib.md5()
    with (gzip.open(file_path, 'rb') if is_compressed else open(file_path, 'rb')) as f:
        for chunk in iter(lambda: f.read(4096), b''):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()

def get_fasta_stats(fasta_file):
    seq_lengths = []
    num_ns = 0
    total_length = 0

    for record in SeqIO.parse(fasta_file, 'fasta'):
        seq_len = len(record)
        seq_lengths.append(seq_len)
        total_length += seq_len
        num_ns += record.seq.lower().count('n')

    sorted_seq_lengths = sorted(seq_lengths, reverse=True)

    def nx(l, x):
        target_len = int(total_length * x)
        cumulative_len = 0
        nx_value = 0
        for length in sorted_seq_lengths:
            cumulative_len += length
            nx_value += 1
            if cumulative_len >= target_len:
                break
        return nx_value

    n50 = nx(sorted_seq_lengths, 0.5)
    n90 = nx(sorted_seq_lengths, 0.9)
    l50 = sorted_seq_lengths[n50 - 1]
    l90 = sorted_seq_lengths[n90 - 1]
    ns_per_100kbp = (num_ns / total_length) * 100000
    num_contigs = len(sorted_seq_lengths)
    largest_contig = sorted_seq_lengths[0]

    return [total_length, n50, n90, l50, l90, ns_per_100kbp, num_contigs, largest_contig]

def main(fasta_file):
    is_compressed = fasta_file.endswith('.gz')
    assembly_id = os.path.splitext(os.path.basename(fasta_file))[0]
    if is_compressed:
        assembly_id = os.path.splitext(assembly_id)[0]

    compressed_md5 = md5(fasta_file, is_compressed) if is_compressed else 'N/A'
    uncompressed_md5 = md5(fasta_file, False)

    if is_compressed:
        with gzip.open(fasta_file, 'rt') as f_in:
            fasta_stats = get_fasta_stats(f_in)
    else:
        with open(fasta_file, 'r') as f_in:
            fasta_stats = get_fasta_stats(f_in)

    output_dir = os.path.dirname(fasta_file)
    output_file = os.path.join(output_dir, assembly_id + "_stats.csv")
    
    header = ['Assembly ID', 'Compressed MD5', 'Uncompressed MD5', 'Total Length', 'N50', 'N90', 'L50', 'L90', '#N per 100kbp', 'Number of Contigs', 'Largest Contig']
    row = [assembly_id, compressed_md5, uncompressed_md5] + fasta_stats

    with open(output_file, 'w', newline='') as f_out:
        csv_writer = csv.writer(f_out)
        csv_writer.writerow(header)
        csv_writer.writerow(row)

    print(f"Statistics saved in {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python fasta_stats.py <input_fasta_file>")
        sys.exit(1)

    fasta_file = sys.argv[1]
    main(fasta_file)
