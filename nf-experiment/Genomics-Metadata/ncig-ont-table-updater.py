# Script Description:
# The script scans a specified directory for ONT genomic data files, 
# extracts detailed metadata from their file paths, and then updates an SQLite database with this information. 
# After updating, it exports the entire database table to an XML file.
# Authors: Kosar Hooshmand

import os
import glob
import sqlite3

# ONT Directory
base_directory_ont = "/g/data/te53/ONT_raw_data/basecalls"
DIR = "/g/data/te53/kh3349/kcgg-data/check"
DB_PATH = os.path.join(DIR, "NCIGSeq.db")

with sqlite3.connect(DB_PATH) as conn:
    cursor = conn.cursor()

    # Process ONT files
    matched_files_ont = sorted(glob.glob(os.path.join(base_directory_ont, '*.fastq.gz')))
    
    for storage_location in matched_files_ont:
        fname = os.path.basename(storage_location)
        parts = fname.split('_')
        
        # Extract metadata from filename
        software_version = "_".join(parts[1:3])  # This captures the 'guppy_x.x.x_hac' as software version

        # Extract sampleID and date_sequenced without considering the repeated sampleID
        sampleID = parts[3]
        date_sequenced = parts[4]
        platform = "ONT"
        
        # Update the NCIG_ids table in the SQLite database
        cursor.execute("""
            UPDATE NCIG_ids 
            SET SoftwareVersion = ?, 
                Sample_ID = ?, 
                DateSequenced = ?, 
                Platform = ?
            WHERE filename1 = ?
        """, (software_version, sampleID, date_sequenced, platform, storage_location))

        if cursor.rowcount == 0:
            print(f"No rows updated for filename1: {storage_location}")
        elif cursor.rowcount > 1:
            print(f"Multiple rows updated for filename1: {storage_location}")

    # Export the entire NCIG_ids table to an XML file
    output_file = os.path.join(DIR, "NCIG_ids_export.xml")

    # Fetch column names for NCIG_ids table
    cursor.execute("PRAGMA table_info(NCIG_ids)")
    columns = [column[1] for column in cursor.fetchall()]

    # Exporting logic
    print(f"Exporting entire table to {output_file}...")
    with open(output_file, "w") as out_f:
        # Starting XML structure
        out_f.write('<?xml version="1.0"?>\n')
        out_f.write('<NCIG_idsTable>\n')
    
        # Writing rows
        for row in cursor.execute("SELECT * FROM NCIG_ids"):
            out_f.write('\t<Row>\n')
            for idx, col in enumerate(columns):
                out_f.write(f'\t\t<{col}>{row[idx]}</{col}>\n')
            out_f.write('\t</Row>\n')
    
        # Closing XML structure
        out_f.write('</NCIG_idsTable>\n')

    print("Exported successfully!")
