import os
import csv
import hashlib
from Bio import SeqIO
import argparse
import gzip


def fasta_to_csv(fasta_file):
    # Extract the assembly ID from the input filename
    assembly_id = os.path.splitext(os.path.basename(fasta_file))[0]

    # Define the output filename based on the assembly ID and input file name
    output_file = f"{os.path.splitext(fasta_file)[0]}_seqtable.csv"

    # Open the output CSV file for writing
    with open(output_file, 'w', newline='') as csvfile:
        # Define the fieldnames for the CSV
        fieldnames = ['Assembly ID', 'Sequence ID', 'Length', 'OrderInFile', 'MD5']
        csv_writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        # Write the header row to the CSV
        csv_writer.writeheader()

        # Determine if the input file is compressed (gzip format)
        is_compressed = False
        if os.path.splitext(fasta_file)[1] == '.gz':
            is_compressed = True

        # Open the input FASTA file for reading
        if is_compressed:
            fasta_handle = gzip.open(fasta_file, 'rt')
        else:
            fasta_handle = open(fasta_file, 'r')

        # Parse the FASTA file
        records = SeqIO.parse(fasta_handle, 'fasta')

        # Iterate over the records in the input FASTA file
        try:
            for idx, record in enumerate(records):
                # Extract the sequence ID and length from the current record
                seq_id = record.id
                seq_length = len(record.seq)

                # Record the order of the current record in the input file
                order_in_file = idx

                # Calculate the MD5 hash of the current sequence
                md5_sum = hashlib.md5(str(record.seq).encode('utf-8')).hexdigest()

                # Create a dictionary representing the current row of the CSV
                row = {
                    'Assembly ID': assembly_id,
                    'Sequence ID': seq_id,
                    'Length': seq_length,
                    'OrderInFile': order_in_file,
                    'MD5': md5_sum
                }

                # Write the current row to the CSV
                csv_writer.writerow(row)
        finally:
            fasta_handle.close()


if __name__ == '__main__':
    # Define a command-line argument parser
    parser = argparse.ArgumentParser(description='Convert a FASTA file to a CSV file with sequence information')

    # Add an argument for the input FASTA file
    parser.add_argument('input_fasta_file', type=str, help='The path to the input FASTA file')

    # Parse the command-line arguments
    args = parser.parse_args()

    # Call the fasta_to_csv function with the input FASTA file
    fasta_to_csv(args.input_fasta_file)
