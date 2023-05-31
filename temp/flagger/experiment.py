from Bio import SeqIO

input_file = "/g/data/xl04/ka6418/species/Tiliqua_Rugosa/rTilRug0.5/assembly/rTilRug0.5.asm.fasta"
output_file = "/g/data/xl04/ka6418/species/Tiliqua_Rugosa/rTilRug0.5/assembly/rTilRug0.5.asm.experimental.fasta"
target_sequence_name = "haplotype1-0000001"
new_sequence_name = "haplotype10"

# Read the FASTA file
sequences = SeqIO.parse(input_file, "fasta")

# Modify the sequence name
modified_sequences = []
for record in sequences:
    if record.id == target_sequence_name:
        record.id = new_sequence_name
        record.description = new_sequence_name
    modified_sequences.append(record)

# Write the modified FASTA file
SeqIO.write(modified_sequences, output_file, "fasta")
