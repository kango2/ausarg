import os
import csv
import hashlib
from Bio import SeqIO
import argparse
import gzip

def fasta_to_csv(fasta_file, output_dir):
    # Extract the assembly ID from the input filename
    assembly_id = os.path.splitext(os.path.basename(fasta_file))[0]

    # Define the output filename based on the assembly ID and input file name
    output_file = os.path.join(output_dir, f"{assembly_id}_seqtable.csv")

    # ... [rest of the function remains unchanged]

if __name__ == '__main__':
    # Define a command-line argument parser
    parser = argparse.ArgumentParser(description='Convert a FASTA file to a CSV file with sequence information')

    # Add an argument for the input FASTA file
    parser.add_argument('-fasta', type=str, required=True, help='The path to the input FASTA file')
    
    # Add an argument for the output directory
    parser.add_argument('-outputdir', type=str, required=True, help='The directory where the output CSV will be placed')

    # Parse the command-line arguments
    args = parser.parse_args()

    # Create output directory if it does not exist
    if not os.path.exists(args.outputdir):
        os.makedirs(args.outputdir)

    # Call the fasta_to_csv function with the input FASTA file and output directory
    fasta_to_csv(args.fasta, args.outputdir)
