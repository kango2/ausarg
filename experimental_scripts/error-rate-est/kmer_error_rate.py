#!/usr/bin/env python3

import sys
import math

def calc_error_rate(histo_file):
    # Read k-mer histogram file
    with open(histo_file, 'r') as f:
        hist_data = [line.strip().split() for line in f.readlines()]

    # Convert histogram data to dictionary
    hist_dict = {}
    for item in hist_data:
        hist_dict[int(item[0])] = int(item[1])

    # Calculate the total number of k-mers
    total_kmers = sum(hist_dict.values())

    # Calculate the number of single-copy k-mers
    single_kmers = hist_dict.get(1, 0)

    # Calculate the number of unique k-mers
    unique_kmers = len(hist_dict)

    # Calculate the estimated error rate
    error_rate = (single_kmers / total_kmers) * math.log(2) / math.log(unique_kmers)

    return error_rate

if __name__ == '__main__':
    # Get input filename from command line
    input_file = sys.argv[1]

    # Calculate error rate
    error_rate = calc_error_rate(input_file)

    # Print error rate to standard output
    print("Estimated error rate: %0.6f" % error_rate)
