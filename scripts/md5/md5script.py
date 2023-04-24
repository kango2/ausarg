import sys
import os
import hashlib
from Bio import SeqIO
import gzip

def md5_checksum(file_path):
    md5_hash = hashlib.md5()
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            md5_hash.update(chunk)
    return md5_hash.hexdigest()

def process_fasta(file_path, is_compressed=False):
    sequences = []

    if is_compressed:
        with gzip.open(file_path, "rt") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                sequences.append((record.id, len(record.seq), hashlib.md5(str(record.seq).encode("utf-8")).hexdigest()))
    else:
        with open(file_path, "rt") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                sequences.append((record.id, len(record.seq), hashlib.md5(str(record.seq).encode("utf-8")).hexdigest()))

    return sequences

def main():
    if len(sys.argv) != 2:
        print("Usage: python script_name.py <FASTA_file>")
        sys.exit(1)

    file_path = sys.argv[1]
    file_name, file_extension = os.path.splitext(file_path)
    is_compressed = file_extension.lower() == ".gz"

    checksum = md5_checksum(file_path)
    if is_compressed:
        uncompressed_checksum = md5_checksum(gzip.open(file_path).name)
    else:
        uncompressed_checksum = None

    sequences = process_fasta(file_path, is_compressed)
    output_file = f"{file_name}_md5.tsv"

    with open(output_file, "w") as f:
        if is_compressed:
            f.write("compressed_md5_checksum\t" + checksum + "\n")
        f.write("uncompressed_md5_checksum\t" + (uncompressed_checksum if is_compressed else checksum) + "\n")
        f.write("\nsequence_name\tlength\tmd5_checksum\n")
        for seq in sequences:
            f.write(f"{seq[0]}\t{seq[1]}\t{seq[2]}\n")

    print(f"Results saved in {output_file}")

if __name__ == "__main__":
    main()
