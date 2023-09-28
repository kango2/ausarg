import argparse
import pandas as pd
from Bio import SeqIO

def count_ns(fasta_file, output_csv):
    data = []

    for record in SeqIO.parse(fasta_file, 'fasta'):
        sequence_name = record.id
        sequence_str = str(record.seq).upper()
        count_n = sequence_str.count('N')
        sequence_length = len(record.seq)
        ratio_n = count_n / sequence_length if sequence_length > 0 else 0
        n_per_100kbp = (count_n / sequence_length * 100000) if sequence_length > 0 else 0
        data.append([sequence_name, count_n, sequence_length, ratio_n, n_per_100kbp])

    df = pd.DataFrame(data, columns=['Sequence Name', 'Count of N', 'Sequence Length', 'Ratio of N', 'N per 100kbp'])
    df = df.sort_values(by='Count of N', ascending=False)
    df.to_csv(output_csv, index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Count Ns in sequences from a FASTA file.')
    parser.add_argument('-i', '--input', required=True, help='Path to the input FASTA file.')
    parser.add_argument('-o', '--output', required=True, help='Path to the output CSV file.')
    args = parser.parse_args()

    count_ns(args.input, args.output)
