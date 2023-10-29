import os
import pysam
import sqlite3
import re

# Constants
CRAM_DIR = "/g/data/te53/chm13_phase3_20220817/"
REFERENCE_FASTA = "/g/data/te53/ontsv/references/chm13_reference_dir/chm13.draft_v1.1.fasta"
DATABASE_PATH = '/g/data/te53/kh3349/kcgg-data/check/NCIGSeq.db'

# Function to determine sex based on X and Y counts
def determine_sex(x_count, y_count):
    ratio = y_count / (x_count + y_count + 1e-5)
    if ratio > 0.2:
        return "Male"
    else:
        return "Female"

# Function to get sex from CRAM file
def get_sex_from_cram(file_path):
    cram_file = pysam.AlignmentFile(file_path, "rc", reference_filename=REFERENCE_FASTA)
    
    x_contig_names = ['X', 'chrX']
    y_contig_names = ['Y', 'chrY']

    x_count = 0
    for contig_name in x_contig_names:
        try:
            x_count += cram_file.count(contig=contig_name)
        except ValueError:
            continue
    
    y_count = 0
    for contig_name in y_contig_names:
        try:
            y_count += cram_file.count(contig=contig_name)
        except ValueError:
            continue

    print(f"X count: {x_count}, Y count: {y_count} for file: {file_path}")  # Diagnostic print
    return determine_sex(x_count, y_count)

# Function to update sex in SQLite database based on extracted info
def update_sex_in_db(sample_id, sex):
    conn = sqlite3.connect(DATABASE_PATH)
    cursor = conn.cursor()
    cursor.execute("UPDATE NCIG_ids SET Sex=? WHERE Sample_ID=?", (sex, sample_id))
    conn.commit()
    conn.close()

# Iterate over all subdirectories in CRAM_DIR for CRAM files
for root, dirs, files in os.walk(CRAM_DIR):
    for file in files:
        if file.endswith('.cram'):
            cram_path = os.path.join(root, file)
            
            # Extract sample_id from the filename
            sample_id_match = re.search(r'FD\d+', file)
            sample_id = sample_id_match.group() if sample_id_match else "N/A"

            try:
                sex = get_sex_from_cram(cram_path)
                print(f"Sample_ID: {sample_id}, Determined Sex: {sex}")  # Diagnostic print
            except Exception as e:
                print(f"Error processing {cram_path}. Error message: {e}. Setting sex to 'null'.")
                sex = 'null'  # Set the sex to null if there's an error

            update_sex_in_db(sample_id, sex)

print("Sex determination and database update completed.")
