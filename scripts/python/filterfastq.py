# FastqHandler created by Hyungtaek Jung
# Filter out based on sequence read quality including butttery-eel basecalled fastq
# Generate two outcomes as passed and failed.fastq

#!/usr/bin/env python3
import argparse
import gzip
import bz2
import logging
import os
import re

# Setup logging
logging.basicConfig(filename='fqqual_filter.log', level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')

def parse_args():
    parser = argparse.ArgumentParser(description='Filter FASTQ files based on sequence read quality or mean_qscore.')
    parser.add_argument('--input-seq', required=True, help='Input FASTQ file and path.')
    parser.add_argument('--filter-qual', type=int, required=True, help='Minimum sequence read quality or mean_qscore to keep (integer).')
    parser.add_argument('--out', required=True, help='Output folder path.')
    parser.add_argument('--t', type=int, default=1, help='Number of CPUs (default: 1).')
    parser.add_argument('--mem', type=int, default=10, help='Amount of memory in GB (default: 10).')
    return parser.parse_args()

def open_file(input_file):
    """Determine the file format and appropriate method to open it."""
    if input_file.lower().endswith('.gz'):
        return gzip.open(input_file, 'rt')
    elif input_file.lower().endswith('.bz2'):
        return bz2.open(input_file, 'rt')
    else:
        return open(input_file, 'r')

def filter_by_quality(input_seq, filter_qual, out_folder):
    """Filter FASTQ sequences based on average read quality or mean_qscore in the header."""
    base_name = os.path.basename(input_seq)  # Get the basename of the input file
    # Strip the file extension and append '_passed' or '_failed'
    base_name_no_ext = os.path.splitext(base_name)[0]
    out_passed = os.path.join(out_folder, f'{base_name_no_ext}_passed.fastq')
    out_failed = os.path.join(out_folder, f'{base_name_no_ext}_failed.fastq')

    with open_file(input_seq) as fq, open(out_passed, 'w') as passed_fq, open(out_failed, 'w') as failed_fq:
        while True:
            header = fq.readline()
            if not header: break  # End of file
            seq = fq.readline()
            plus = fq.readline()
            qual = fq.readline()
            mean_qscore_match = re.search(r'mean_qscore=(\d+)', header)
            if mean_qscore_match:
                avg_qual = int(mean_qscore_match.group(1))
            else:
                qual_scores = [ord(char) - 33 for char in qual.strip()]
                avg_qual = sum(qual_scores) / len(qual_scores)
            if avg_qual >= filter_qual:
                passed_fq.write(header + seq + plus + qual)
            else:
                failed_fq.write(header + seq + plus + qual)

def main():
    args = parse_args()

    if args.t > 1:
        logging.warning("Multiprocessing is noted but not implemented due to script design.")

    if not os.path.exists(args.out):
        os.makedirs(args.out)

    filter_by_quality(args.input_seq, args.filter_qual, args.out)

    logging.info(f"Quality filtering complete. Outputs written to {args.out}")

if __name__ == '__main__':
    main()

