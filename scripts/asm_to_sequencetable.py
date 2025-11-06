import os
import csv
import hashlib
from Bio import SeqIO
import argparse
import gzip
from concurrent.futures import ProcessPoolExecutor

def calc_avg_gc(seq, window=10000):
    if len(seq) == 0:
        return 0.0
    gc_total = 0
    window_count = 0
    for i in range(0, len(seq) - window + 1, window):
        subseq = seq[i:i+window].upper()
        gc = subseq.count('G') + subseq.count('C')
        gc_total += (gc / window) * 100
        window_count += 1
    tail_start = (len(seq) // window) * window
    if tail_start < len(seq):
        tail = seq[tail_start:].upper()
        gc = tail.count('G') + tail.count('C')
        gc_total += (gc / len(tail)) * 100
        window_count += 1
    return round(gc_total / window_count, 2) if window_count else 0.0

def process_record(args):
    record, assembly_id = args
    seq_id = record.id
    seq = str(record.seq)
    seq_length = len(seq)
    md5_sum = hashlib.md5(seq.encode('utf-8')).hexdigest()
    avg_gc = calc_avg_gc(seq)
    return {
        'asmid': assembly_id,
        'seqid': seq_id,
        'len': seq_length,
        'md5': md5_sum,
        'gc': avg_gc
    }

def fasta_to_csv_parallel(fasta_file, output_dir, threads):
    assembly_id = os.path.splitext(os.path.basename(fasta_file))[0]
    output_file = os.path.join(output_dir, f"{assembly_id}_seqtable.csv")

    is_compressed = fasta_file.endswith('.gz')
    fasta_handle = gzip.open(fasta_file, 'rt') if is_compressed else open(fasta_file, 'r')
    records = list(SeqIO.parse(fasta_handle, 'fasta'))
    fasta_handle.close()

    args_list = [(record, assembly_id) for record in records]

    with ProcessPoolExecutor(max_workers=threads) as executor:
        results = list(executor.map(process_record, args_list))

    with open(output_file, 'w', newline='') as csvfile:
        fieldnames = ['asmid', 'seqid', 'len', 'md5', 'gc']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for row in results:
            writer.writerow(row)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate a sequence table with MD5 and GC% from FASTA (parallelized)')
    parser.add_argument('-fasta', required=True, help='Input FASTA file (.fasta or .fasta.gz)')
    parser.add_argument('-outputdir', required=True, help='Directory to write output CSV')
    parser.add_argument('-p', '--processes', type=int, default=4, help='Number of parallel processes')
    args = parser.parse_args()

    os.makedirs(args.outputdir, exist_ok=True)
    fasta_to_csv_parallel(args.fasta, args.outputdir, args.processes)
