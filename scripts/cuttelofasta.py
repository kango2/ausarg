import argparse
from Bio import SeqIO
from Bio.Seq import Seq

def modify_fasta(input_file, output_file):
    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        for record in SeqIO.parse(infile, "fasta"):
            sequence_length = len(record.seq)
            if sequence_length > 80000:
                # Retain first 30,000 and last 30,000 base pairs, replace the middle region with Ns
                retained_start = record.seq[:30000]
                retained_end = record.seq[-30000:]
                middle_ns = 'N' * (sequence_length - 60000)
                new_seq = retained_start + middle_ns + retained_end
            else:
                # Retain the entire sequence if it is less than or equal to 80,000 base pairs
                new_seq = record.seq
            # Write the modified sequence to the output file
            record.seq = Seq(new_seq)
            SeqIO.write(record, outfile, "fasta")

def main():
    parser = argparse.ArgumentParser(description="Modify a FASTA file to retain only the first 30,000 and last 30,000 base pairs of each sequence if the sequence is longer than 80,000 base pairs, replacing the middle region with Ns.")
    parser.add_argument('-input', required=True, help="Path to the input FASTA file")
    parser.add_argument('-output', required=True, help="Path to the output FASTA file")

    args = parser.parse_args()
    
    modify_fasta(args.input, args.output)

if __name__ == "__main__":
    main()
