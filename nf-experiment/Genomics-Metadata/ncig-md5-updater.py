# Script Description:
# This script checks the integrity of genomic data files by verifying the presence of their MD5 checksums. 
# It searches specific directories for genomic files and their associated .md5 files. 
# If a matching .md5 file is found for a genomic file, the script updates an SQLite database, marking that file's integrity as verified. 
# After processing, the updated database information is exported to a .tsv file.
# Authors: Kosar Hooshmand

import sqlite3
import os

# Specify the directory
DIR = "/g/data/te53/kh3349/kcgg-data"
DB_PATH = os.path.join(DIR, "NCIGSeq.db")

# Connect to the SQLite database
conn = sqlite3.connect(DB_PATH)
cursor = conn.cursor()

# Search for subdirectories in the specified path
directories = [d for d in os.walk("/g/data/te53/kccg_illumina_hiseqx10_wgs/") if os.path.isdir(d[0])]

# Loop over directories
for dirpath, _, _ in directories:
    # Find all *R1.fastq.gz and *R.fastq.gz files in the current directory
    for root, _, files in os.walk(dirpath):
        for file in files:
            if file.endswith("_R1.fastq.gz") or file.endswith("_R1_*.fastq.gz"):
                path = os.path.join(root, file)
                # Check if its .md5 counterpart exists
                if os.path.exists(path + ".md5"):
                    # Update the database using the full path
                    cursor.execute("UPDATE NCIG_ids SET Md5Sum = 'exist' WHERE filename1 = ?", (path,))
                    conn.commit()

# Check if the file transfer_from_ontsv_md5sum.txt exists
txt_path = "/g/data/te53/ONT_raw_data/basecalls/transfer_from_ontsv_md5sum.txt"
if os.path.exists(txt_path):
    with open(txt_path, "r") as f:
        for line in f:
            md5, file = line.strip().split()  # Assuming space-separated values
            path = os.path.join("/g/data/te53/ONT_raw_data/basecalls/", file)
            # Check if the file exists
            if os.path.exists(path):
                # Update the NCIG_ids table in NCIGSeq.db
                cursor.execute("UPDATE NCIG_ids SET Md5Sum = 'exist' WHERE filename1 = ?", (path,))
                conn.commit()

# Export the updated table to a .tsv file
output_file = "NCIG_ids_md5sum.tsv"

# Fetching column names for NCIG_ids table
cursor.execute("PRAGMA table_info(NCIG_ids)")
columns = [column[1] for column in cursor.fetchall()]

# Exporting logic
print(f"Exporting updated table to {output_file}...")
with open(output_file, "w") as out_f:
    # Writing column headers
    out_f.write("\t".join(columns) + "\n")
    
    # Writing rows
    for row in cursor.execute("SELECT * FROM NCIG_ids"):
        if row[columns.index('Md5Sum')] is not None:
            out_f.write("\t".join(map(str, row)) + "\n")

print("Done!")
conn.close()
