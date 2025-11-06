import argparse
import os
from Bio import SeqIO
from concurrent.futures import ProcessPoolExecutor

def find_n_regions(seq_record):
    start = None
    regions = []
    for idx, base in enumerate(seq_record.seq):
        if base.upper() == 'N':
            if start is None:
                start = idx
        else:
            if start is not None:
                # BED format: 0-based start, 1-based end
                regions.append((seq_record.id, start, idx, idx - start))
                start = None
    if start is not None:
        regions.append((seq_record.id, start, len(seq_record), len(seq_record) - start))
    
    return regions  # return empty list if no N-regions

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find N-gap regions in a FASTA and output as BED.")
    parser.add_argument("-i", "--input", required=True, help="Path to the FASTA file.")
    parser.add_argument("-o", "--output", required=True, help="Path to the output directory.")
    parser.add_argument("-p", "--processes", type=int, default=4, help="Number of parallel processes.")
    
    args = parser.parse_args()
    
    output_file_path = os.path.join(
        args.output,
        os.path.basename(args.input).rsplit('.', 1)[0] + '_gaps.bed'
    )

    with open(args.input, 'r') as f, open(output_file_path, 'w') as out_f:
        sequences = list(SeqIO.parse(f, "fasta"))
        with ProcessPoolExecutor(max_workers=args.processes) as executor:
            results = executor.map(find_n_regions, sequences)
            
            for seq_regions in results:
                for chrom, start, end, width in seq_regions:
                    out_f.write(f"{chrom}\t{start}\t{end}\t{width}\n")