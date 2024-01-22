import argparse
from Bio import SeqIO

def calculate_gc_content(seq, window_size, overlap):
    gc_content = []
    for i in range(0, len(seq) - window_size + 1, window_size - overlap):
        window = seq[i:i + window_size].upper()  # Convert to uppercase
        gc_count = window.count('G') + window.count('C')
        gc_content.append((i, i + window_size, gc_count / window_size))
    return gc_content

def main():
    parser = argparse.ArgumentParser(description='Calculate GC content in a sliding window for a FASTA file and output in BED format.')
    parser.add_argument('input_file', type=str, help='Input FASTA file')
    parser.add_argument('output_file', type=str, help='Output file to save the results')
    parser.add_argument('window_size', type=int, help='Size of the sliding window')
    parser.add_argument('overlap', type=int, help='Overlap size between windows')

    args = parser.parse_args()

    with open(args.input_file, "r") as fasta_file, open(args.output_file, "w") as output_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            gc_content = calculate_gc_content(str(record.seq), args.window_size, args.overlap)
            for start, end, gc in gc_content:
                output_file.write(f"{record.id}\t{start}\t{end}\t{gc:.2f}\n")

if __name__ == "__main__":
    main()
