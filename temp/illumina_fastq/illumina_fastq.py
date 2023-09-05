import csv
import gzip
from collections import defaultdict
import hashlib
import argparse
import os
import json

def parse_fastq(fastq_file):
    """Parse the FASTQ file and return a list of sequences, quality scores, and MD5 checksums."""
    sequences = []
    quality_scores = []
    hash_md5_compressed = hashlib.md5()
    hash_md5_uncompressed = hashlib.md5()

    with open(fastq_file, 'rb') as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5_compressed.update(chunk)

    with gzip.open(fastq_file, 'rt') as f:
        for i, line in enumerate(f):
            # Sequence lines are the second line in every set of 4 lines
            if i % 4 == 1:
                sequences.append(line.strip())
            # Quality scores are the fourth line in every set of 4 lines
            elif i % 4 == 3:
                quality_scores.append(line.strip())

            # Update the MD5 hash for the uncompressed content
            hash_md5_uncompressed.update(line.encode())

    md5_compressed = hash_md5_compressed.hexdigest()
    md5_uncompressed = hash_md5_uncompressed.hexdigest()

    return sequences, quality_scores, md5_compressed, md5_uncompressed

def count_reads(sequences):
    """Count the number of reads from the sequences."""
    return len(sequences)

def count_total_bases(sequences):
    """Count the total number of bases from the sequences."""
    return sum(len(seq) for seq in sequences)

def average_read_length(num_reads, total_bases):
    """Calculate the average read length."""
    return total_bases / num_reads if num_reads > 0 else 0

def phred_to_quality(phred_string):
    """Convert Phred encoded quality string to a list of numerical quality values."""
    return [ord(char) - 33 for char in phred_string]

def average_per_base_quality(quality_scores):
    """Calculate average per base quality values across all reads."""
    total_scores = [0] * max(len(qs) for qs in quality_scores)
    count_scores = [0] * max(len(qs) for qs in quality_scores)  # Count of scores for each position
    
    for qs in quality_scores:
        scores = phred_to_quality(qs)
        for i, score in enumerate(scores):
            total_scores[i] += score
            count_scores[i] += 1  # Increment count for that position

    # Calculate average by dividing total score by count for each position
    return [round(total_scores[i] / count_scores[i]) if count_scores[i] != 0 else 0 for i in range(len(total_scores))]

def average_nucleotide_frequencies(sequences):
    """Calculate average per base nucleotide frequencies (A, C, G, T) across all reads."""
    bases = ['A', 'C', 'G', 'T']
    total_bases = [defaultdict(int) for _ in range(max(len(seq) for seq in sequences))]
    total_reads = len(sequences)
    
    for seq in sequences:
        for i, base in enumerate(seq):
            total_bases[i][base] += 1

    average_frequencies = []
    for base_counts in total_bases:
        avg_freq = [base_counts.get(base, 0) for base in bases]
        average_frequencies.append(avg_freq)

    return average_frequencies

def overall_nucleotide_content(sequences):
    """Calculate overall nucleotide content (A, C, G, T) for the entire FASTQ file."""
    bases = ['A', 'C', 'G', 'T']
    total_counts = defaultdict(int)
    total_length = sum(len(seq) for seq in sequences)
    
    for seq in sequences:
        for base in seq:
            total_counts[base] += 1

    return [total_counts.get(base, 0) for base in bases]

def average_qv_per_read(quality_scores):
    """Calculate histogram of average quality values binned in intervals of 1."""
    
    # Initialize the histogram with 101 bins (0-100)
    histogram = [0] * 101
    
    for qs in quality_scores:
        scores = phred_to_quality(qs)
        
        # Calculate the average quality for this read
        avg_quality = sum(scores) / len(scores) if scores else 0
        
        # Determine the bin for this average quality
        bin_index = min(int(avg_quality), 100)  # max bin index is 100
        
        # Increment the count for the determined bin
        histogram[bin_index] += 1

    return histogram

def calculate_gc_content(sequences):
    """Calculate the GC content percentage for the entire FASTQ file."""
    total_length = sum(len(seq) for seq in sequences)
    gc_count = sum(seq.count('G') + seq.count('C') for seq in sequences)
    gc_content = (gc_count / total_length) * 100 if total_length > 0 else 0
    return round(gc_content, 2)

def get_file_size_gb(file_path):
    """Return the file size in gigabytes."""
    size_in_bytes = os.path.getsize(file_path)
    size_in_gb = size_in_bytes / (1024 ** 3)
    return round(size_in_gb, 2)

def main():
    parser = argparse.ArgumentParser(description='Calculate metrics for Illumina FASTQ files (R1 and R2).')
    parser.add_argument('-R1', type=str, required=True, help='Path to the R1 fastq.gz file.')
    parser.add_argument('-R2', type=str, required=True, help='Path to the R2 fastq.gz file.')
    parser.add_argument('-o', '--output', type=str, default="statistics.csv", help='Path to the output CSV file.')
    args = parser.parse_args()

    metrics = []

    for fastq_file in [args.R1, args.R2]:
        sequences, quality_scores, md5_compressed, md5_uncompressed = parse_fastq(fastq_file)
        num_reads = count_reads(sequences)
        total_bases = count_total_bases(sequences)
        avg_read_length_val = average_read_length(num_reads, total_bases)
        avg_quality_values = average_per_base_quality(quality_scores)
        nucleotide_freq = average_nucleotide_frequencies(sequences)
        overall_content = overall_nucleotide_content(sequences)
        avg_qv_reads = average_qv_per_read(quality_scores)
        gc_content = calculate_gc_content(sequences)
        file_size = get_file_size_gb(fastq_file)

        # Format nucleotide frequencies
        nucleotide_freq_str = [":".join(str(int(freq)) if float(freq).is_integer() else f"{freq:.2f}" for freq in freq_list) for freq_list in nucleotide_freq]
        overall_content_str = ":".join(f"{freq:.2f}" for freq in overall_content)

        metrics.append({
            'num_reads': num_reads,
            'total_bases': total_bases,
            'avg_read_length': avg_read_length_val,
            'avg_quality_values': avg_quality_values,
            'nucleotide_freq': nucleotide_freq_str,
            'overall_content': overall_content_str,
            'avg_qv_reads': avg_qv_reads,
            'md5_compressed': md5_compressed,
            'md5_uncompressed': md5_uncompressed,
            'gc_content': gc_content,
            'file_size_gb': file_size
        })

    output_csv = args.output
    with open(output_csv, 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        headers = ['Number of reads', 'Number of bases', 'Mean read length', 
                   'Mean QV at read position', 'Nucleotide count at read position', 
                   'Nucleotide content ', 'Mean QV per read',
                   'MD5 (.fastq.gz)', 'MD5 (.fastq)', 'GC', 'File size in GB']
        csvwriter.writerow(headers)

        combined_metrics = []
        keys = ['num_reads', 'total_bases', 'avg_read_length', 'avg_quality_values',
                'nucleotide_freq', 'overall_content', 'avg_qv_reads', 'md5_compressed', 
                'md5_uncompressed', 'gc_content', 'file_size_gb']  # Added gc_content and file_size_gb


        for key in keys:
            if isinstance(metrics[0][key], list):
                combined_val = ",".join(str(val) for val in metrics[0][key]) + ";" + ",".join(str(val) for val in metrics[1][key])
            else:
                combined_val = f"{metrics[0][key]};{metrics[1][key]}"
            combined_metrics.append(combined_val)
        
        csvwriter.writerow(combined_metrics)
    
    metrics_dict = dict(zip(headers, combined_metrics))
    
    # Write to JSON
    output_json = args.output.replace('.csv', '.json')
    with open(output_json, 'w') as jsonfile:
        json.dump({"metrics": metrics_dict}, jsonfile, indent=4)
    
    print(f"Statistics saved to {output_csv} and {output_json}")

if __name__ == "__main__":
    main()