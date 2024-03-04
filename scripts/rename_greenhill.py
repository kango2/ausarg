# Path to your original FASTA file
original_fasta_file_path = '/g/data/xl04/ka6418/bassiana/publication/asm_fasta/BASDU_HifiASM_GreenHill_SUP_renamed.fasta'

# Path for the new FASTA file with sequences renamed in descending order of their length
new_fasta_file_path = '/g/data/xl04/ka6418/bassiana/publication/asm_fasta/BASDU_HifiASM_GreenHill_SUP_renamed_sequences.fasta'

# Function to read and store sequences from the original FASTA file
def read_fasta(file_path):
    sequences = []
    with open(file_path, 'r') as file:
        sequence_name = ''
        sequence_data = ''
        for line in file:
            if line.startswith('>'):
                if sequence_data:
                    sequences.append((sequence_name, sequence_data, len(sequence_data)))
                    sequence_data = ''
                sequence_name = line.strip()
            else:
                sequence_data += line.strip()
        if sequence_data:
            sequences.append((sequence_name, sequence_data, len(sequence_data)))
    return sequences

# Read sequences from the original FASTA file
sequences = read_fasta(original_fasta_file_path)

# Sort the sequences by length in descending order
sorted_sequences = sorted(sequences, key=lambda x: x[2], reverse=True)

# Write the renamed and sorted sequences to the new FASTA file
with open(new_fasta_file_path, 'w') as new_file:
    sequence_number = 1
    for _, sequence_data, _ in sorted_sequences:
        new_file.write(f'>seq{sequence_number}\n')
        new_file.write(f'{sequence_data}\n')
        sequence_number += 1

print(f'Renamed and reordered sequences FASTA file saved to: {new_fasta_file_path}')
