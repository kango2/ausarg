import csv
import gzip
from collections import defaultdict
import hashlib
import argparse
import os
import json
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor

def phred_to_quality(phred_string):
    """Convert Phred encoded quality string to a list of numerical quality values."""
    return [ord(char) - 33 for char in phred_string]

def get_file_size_gb(file_path):
    """Return the file size in gigabytes."""
    size_in_bytes = os.path.getsize(file_path)
    size_in_gb = size_in_bytes / (1024 ** 2)  # Corrected to return in GB
    return round(size_in_gb, 2)

def calculate_md5_compressed(fastq_file):
    hash_md5_compressed = hashlib.md5()
    with open(fastq_file, 'rb') as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5_compressed.update(chunk)
    return hash_md5_compressed.hexdigest()

def calculate_md5_uncompressed(fastq_file):
    hash_md5_uncompressed = hashlib.md5()
    with gzip.open(fastq_file, 'rt') as f:
        for line in f:
            hash_md5_uncompressed.update(line.encode())
    return hash_md5_uncompressed.hexdigest()

def parse_fastq(fastq_file):
    """Parse the FASTQ file and return metrics."""
    num_reads = 0
    total_bases = 0
    quality_scores_sum = defaultdict(int)
    quality_scores_count = defaultdict(int)
    total_bases_count = defaultdict(int)
    avg_qv_reads = [0] * 101
    total_gc_count = 0
    nucleotide_freq_data = defaultdict(lambda: defaultdict(int))
    overall_content_data = defaultdict(int)
    bases = set()
    with gzip.open(fastq_file, 'rt') as f:
        for i, line in enumerate(f):
            # Sequence lines are the second line in every set of 4 lines
            if i % 4 == 1:
                seq = line.strip()
                num_reads += 1
                total_bases += len(seq)
                total_gc_count += seq.count('G') + seq.count('C')
                for idx, base in enumerate(seq):
                    total_bases_count[idx] += 1
                    nucleotide_freq_data[idx][base] += 1
                    overall_content_data[base] += 1
                    bases.update(base)
            elif i % 4 == 3:
                qs = line.strip()
                scores = phred_to_quality(qs)
                avg_quality = sum(scores) / len(scores) if scores else 0
                bin_index = min(int(avg_quality), 100)
                avg_qv_reads[bin_index] += 1
                for idx, score in enumerate(scores):
                    quality_scores_sum[idx] += score
                    quality_scores_count[idx] += 1

    avg_quality_values = [round(quality_scores_sum[i] / quality_scores_count[i]) if quality_scores_count[i] != 0 else 0 for i in range(len(quality_scores_sum))]
    gc_content = round((total_gc_count / total_bases) * 100, 2) if total_bases else 0
    nucleotide_freq = [":".join(str(nucleotide_freq_data[idx].get(base, 0)) for base in bases) for idx in nucleotide_freq_data]
    overall_content = ":".join(str(overall_content_data.get(base, 0)) for base in bases)

    with ThreadPoolExecutor() as executor:
        future_md5_compressed = executor.submit(calculate_md5_compressed, fastq_file)
        future_md5_uncompressed = executor.submit(calculate_md5_uncompressed, fastq_file)
        md5_compressed = future_md5_compressed.result()
        md5_uncompressed = future_md5_uncompressed.result()

    return {
        'num_reads': num_reads,
        'total_bases': total_bases,
        'avg_read_length': total_bases / num_reads if num_reads > 0 else 0,
        'avg_quality_values': avg_quality_values,
        'gc_content': gc_content,
        'md5_compressed': md5_compressed,
        'md5_uncompressed': md5_uncompressed,
        'avg_qv_reads': avg_qv_reads,
        'file_size_gb': get_file_size_gb(fastq_file),
        'path': os.path.abspath(fastq_file),
        'nucleotide_freq': nucleotide_freq,
        'overall_content': overall_content
    }

def main():
    parser = argparse.ArgumentParser(description='Calculate metrics for Illumina FASTQ files (R1 and R2).')
    parser.add_argument('-R1', type=str, required=True, help='Path to the R1 fastq.gz file.')
    parser.add_argument('-R2', type=str, help='Path to the R2 fastq.gz file. Optional.')
    parser.add_argument('-o', '--output', type=str, default="statistics", help='Base name for the output file (without extension).')
    parser.add_argument('--format', type=str, choices=["csv", "json", "both"], default="csv", help='Output format: csv, json, or both.')
    parser.add_argument('--cores', type=int, default=os.cpu_count(), help='Number of cores to use. Defaults to all available cores.')
    args = parser.parse_args()

    metrics = []
    fastq_files = [args.R1]
    if args.R2:
        fastq_files.append(args.R2)

    num_cores_per_file = args.cores // len(fastq_files)

    with ProcessPoolExecutor(max_workers=num_cores_per_file) as executor:
        for result in executor.map(parse_fastq, fastq_files):
            metrics.append(result)

    headers = [
        'File_path','Number_of_reads', 'Number_of_bases', 'Mean_read_length',
        'Mean_QV_at_read_position', 'Nucleotide_count_at_read_position',
        'Nucleotide_content', 'Mean_QV_per_read',
        'MD5_zipped', 'MD5_text', 'GC', 'File_size_in_MB'
    ]

    combined_metrics = []
    keys = [
        'path','num_reads', 'total_bases', 'avg_read_length', 'avg_quality_values',
        'nucleotide_freq', 'overall_content', 'avg_qv_reads', 'md5_compressed',
        'md5_uncompressed', 'gc_content', 'file_size_gb'
    ]

    for key in keys:
        if isinstance(metrics[0][key], list):
            combined_val = ",".join(str(val) for val in metrics[0][key])
            if len(metrics) > 1:
                combined_val += ";" + ",".join(str(val) for val in metrics[1][key])
        else:
            combined_val = f"{metrics[0][key]}"
            if len(metrics) > 1:
                combined_val += f";{metrics[1][key]}"
        combined_metrics.append(combined_val)

    metrics_dict = dict(zip(headers, combined_metrics))

    if "both" in args.format or "csv" in args.format:
        with open(args.output + '.csv', 'w', newline='') as csvfile:
            csvwriter = csv.writer(csvfile)
            csvwriter.writerow(headers)
            csvwriter.writerow(combined_metrics)

    if "both" in args.format or "json" in args.format:
        output_json = args.output.replace('.csv', '.json')
        with open(output_json + '.json', 'w') as jsonfile:
            json.dump({"metrics": metrics_dict}, jsonfile, indent=4)

if __name__ == "__main__":
    main()
