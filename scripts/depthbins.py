import csv
import argparse
from collections import defaultdict
import os

def parse_arguments():
    parser = argparse.ArgumentParser(description='Calculate read depth frequencies binned by range.')
    parser.add_argument('-input', type=str, required=True, help='Path to the input CSV file.')
    parser.add_argument('-output', type=str, required=False, help='Path to the output file. Optional.')
    return parser.parse_args()

def find_bin(read_depth, bins):
    """Determine the appropriate bin for a given read depth, capping at 200."""
    if read_depth > 200:
        return "190-200"
    for i in bins[:-1]:
        if i <= read_depth < i + 10:
            return f"{i}-{i+10}"
    return "190-200"

def main():
    args = parse_arguments()
    input_file = args.input
    output_file = args.output if args.output else os.path.splitext(input_file)[0] + "_depth_bins.txt"

    bins = range(0, 201, 10)
    bin_labels = [f"{i}-{i+10}" for i in bins[:-1]]
    read_depth_frequency = defaultdict(int)

    with open(input_file, 'r') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            read_depth = float(row[3])
            appropriate_bin = find_bin(read_depth, bins)
            read_depth_frequency[appropriate_bin] += 1

    for label in bin_labels:
        if label not in read_depth_frequency:
            read_depth_frequency[label] = 0

    sorted_frequencies = sorted(read_depth_frequency.items(), key=lambda x: int(x[0].split("-")[0]))

    with open(output_file, 'w') as f:
        for bin_range, frequency in sorted_frequencies:
            f.write(f"{bin_range}\t{frequency}\n")

    print(f"Output saved to {output_file}")

if __name__ == "__main__":
    main()
