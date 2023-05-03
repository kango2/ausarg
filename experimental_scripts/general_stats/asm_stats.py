#!/usr/bin/python

"""
This script calculates all the basic length metrics for the given input genome assembly.

Usage: python3/Python-3.6.4/python assembly_stats.py <genome assembly file (fasta format)> <output file name> <estimated genome size (in Mb)> 

"""

from Bio import SeqIO
import sys
import statistics
import numpy as np

inputfile = sys.argv[1]
outputfile = sys.argv[2]
estimated_genome_size = float(sys.argv[3])

# Initialize a dictionary to store the variable names and values
assembly_stats = {}

records = list(SeqIO.parse(inputfile, "fasta"))

number_of_scaffolds = len(records)
assembly_stats["Number_of_scaffolds"] = number_of_scaffolds

len_seq = [len(rec) for rec in records]

total_size_scaffolds = sum(len_seq)
assembly_stats["Total_size_scaffolds"] = total_size_scaffolds

total_scaffold_length_percentage_genome_size = ((total_size_scaffolds/(estimated_genome_size*1000000))*100)
assembly_stats["Total_scaffold_length_percentage_genome_size"] = total_scaffold_length_percentage_genome_size



# Write the dictionary to a TSV file
with open(outputfile, "w") as output:
    output.write("Metric\tValue\n")
    for key, value in assembly_stats.items():
        output.write("{}\t{}\n".format(key, value))
