import os
import pysam
import sqlite3

# Function to determine sex based on X and Y counts
def determine_sex(x_count, y_count):
    ratio = y_count / (x_count + y_count + 1e-5)
    if ratio > 0.2:
        return "Male"
    else:
        return "Female"

# Function to get sex from BAM file
def get_sex_from_bam(file_path):
    bam_file = pysam.AlignmentFile(file_path, "rb")
    
    x_contig_names = ['X', 'chrX']
    y_contig_names = ['Y', 'chrY']

    x_count = 0
    for contig_name in x_contig_names:
        try:
            x_count += bam_file.count(contig=contig_name)
        except ValueError:
            continue
    
    y_count = 0
    for contig_name in y_contig_names:
        try:
            y_count += bam_file.count(contig=contig_name)
        except ValueError:
            continue

    print(f"X count: {x_count}, Y count: {y_count} for file: {file_path}")  # Diagnostic print
    return determine_sex(x_count, y_count)

# Function to update sex in SQLite database
def update_sex_in_db(sample_id, sex):
    DATABASE_PATH = '/g/data/te53/kh3349/kcgg-data/check/NCIGSeq.db'  # <-- Replace with the actual path to your database
    conn = sqlite3.connect(DATABASE_PATH)
    cursor = conn.cursor()
    cursor.execute("UPDATE NCIG_ids SET Sex=? WHERE Sample_ID=?", (sex, sample_id))
    conn.commit()
    conn.close()

# Specify BAM directory
BAM_DIR = "/g/data/te53/ontsv/data/alignments/minimap2/v2.22-r1101"

# Iterate over all subdirectories in BAM_DIR
for root, dirs, files in os.walk(BAM_DIR):
    for file in files:
        if file.endswith('_pass.bam'):
            print(f"Found BAM file: {os.path.join(root, file)}")  # Diagnostic print
            bam_path = os.path.join(root, file)
            sample_id = os.path.basename(bam_path).replace('_pass.bam', '')
            sex = get_sex_from_bam(bam_path)
            print(f"Sample: {sample_id}, Determined Sex: {sex}")  # Diagnostic print
            update_sex_in_db(sample_id, sex)

print("Sex determination and database update completed.")
