#!/usr/bin/env python3

import sys

def main():
    if len(sys.argv) != 3:
        print("Usage: python process_bedgraph.py input_file.bedgraph bin_size")
        sys.exit(1)

    filename = sys.argv[1]
    bin_size = int(sys.argv[2])

    with open(filename) as bedgraph:
        for line in bedgraph:
            fields = line.strip().split('\t')
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            value = fields[3]
            if end - start > bin_size:
                num_bins = (end - start) // bin_size
                if (end - start) % bin_size != 0:
                    num_bins += 1
                for i in range(num_bins):
                    bin_start = start + i * bin_size
                    bin_end = min(end, start + (i+1) * bin_size)
                    print(f"{chrom}\t{bin_start}\t{bin_end}\t{value}")
            else:
                print(line.strip())

if __name__ == '__main__':
    main()