import hashlib
from Bio import SeqIO
import gzip

def md5_checksum(data):
    md5 = hashlib.md5()
    md5.update(data)
    return md5.hexdigest()

input_file = 'rTilRug0.5.asm.hp1.1.fasta'
output_file = 'output.md5.txt'

# Calculate the md5 checksum for the zipped file
with open(input_file, 'rb') as input:
    zipped_checksum = md5_checksum(input.read())

# Calculate the md5 checksum for the unzipped file
if input_file.endswith('.gz'):
    with gzip.open(input_file, 'rb') as input, open(output_file, 'w') as output:
        unzipped_checksum = md5_checksum(input.read())
        output.write(f'Zipped checksum: {zipped_checksum}\n')
        output.write(f'Unzipped checksum: {unzipped_checksum}\n')
else:
    with open(input_file, 'rb') as input, open(output_file, 'w') as output:
        unzipped_checksum = md5_checksum(input.read())
        output.write(f'Unzipped checksum: {unzipped_checksum}\n')

# Calculate the md5 checksum for each sequence in the FASTA file
with open(input_file, 'r') as fasta_in, open(output_file, 'a') as output:
    for record in SeqIO.parse(fasta_in, 'fasta'):
        checksum = md5_checksum(str(record.seq).encode('utf-8'))
        seq_length = len(record.seq)
        output.write(f'{record.id}\t{seq_length}\t{checksum}\n')
