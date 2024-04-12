import argparse
from collections import defaultdict
import os

def parse_arguments():
    parser = argparse.ArgumentParser(description='Calculate read depth frequencies for each value up to 200.')
    parser.add_argument('-input', type=str, required=True, help='Path to the input BED file.')
    parser.add_argument('-output', type=str, required=False, help='Path to the output file. Optional.')
    return parser.parse_args()

def main():
    args = parse_arguments()
    input_file = args.input
    output_file = args.output if args.output else os.path.splitext(input_file)[0] + "_depth_frequencies.txt"

    read_depth_frequency = defaultdict(int)

    with open(input_file, 'r') as bedfile:
        for line in bedfile:
            columns = line.strip().split()
            try:
                read_depth = round(float(columns[3]))  # Assuming the fourth column is the read depth
                if read_depth <= 200:
                    read_depth_frequency[read_depth] += 1
            except ValueError:
                # Handle the case where conversion to float fails
                print(f"Warning: Skipping invalid read depth value '{columns[3]}'")

    # Ensure all values up to 200 are accounted for in the output
    for depth in range(201):  # Include 0 through 200
        if depth not in read_depth_frequency:
            read_depth_frequency[depth] = 0

    # Sort the frequencies by read depth
    sorted_frequencies = sorted(read_depth_frequency.items())

    with open(output_file, 'w') as f:
        for read_depth, frequency in sorted_frequencies:
            f.write(f"{read_depth}\t{frequency}\n")

    print(f"Output saved to {output_file}")

if __name__ == "__main__":
    main()
