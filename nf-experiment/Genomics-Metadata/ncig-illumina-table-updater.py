# Script Description:
# The script scans a specified directory for genomic data files, extracts detailed metadata from their file paths, and then updates an SQLite database with 
# this information. 
# After updating, it exports the entire database table to a .tsv file.
# Authors: Kosar Hooshmand

import os
import glob
import sqlite3

base_directory = "/g/data/te53/kccg_illumina_hiseqx10_wgs/"
DIR = "/g/data/te53/kh3349/kcgg-data/check"
DB_PATH = os.path.join(DIR, "NCIGSeq.db")

with sqlite3.connect(DB_PATH) as conn:
    cursor = conn.cursor()
    matched_files = sorted(glob.glob(os.path.join(base_directory, '**', '*R1.fastq.gz'), recursive=True))

    for storage_location in matched_files:
        fname = os.path.basename(storage_location)
        parts = fname.split('_')
        
        # Determine if the current file is from the "RenalFailureStudy"
        is_renal_failure_study = "RenalFailureStudy" in storage_location

        if is_renal_failure_study:
            flowcellid = parts[0]
            lane = parts[1]
            date_sequenced = parts[2]
            sampleID = parts[3]
            species = parts[4]
            barcode = None  # Directly assigning None since it's missing
            runIdentifier = parts[6]
            technician = parts[7]
            batch = parts[8]
            run = parts[9]
            date_collected = None  # No date_collected for RenalFailureStudy
        else:
            flowcellid = parts[0]
            lane = parts[1]
            date_sequenced = parts[2]
            sampleID = parts[3]
            species = parts[4]
            
            # Locate "R" and assign values after that
            runIdentifier_index = parts.index("R")
            runIdentifier = "R"
            
            # Assign date_collected only if not from "RenalFailureStudy"
            date_collected = parts[runIdentifier_index + 1]
            
            # Check if there is a Barcode
            if "-" in parts[runIdentifier_index - 1]:
                barcode = parts[runIdentifier_index - 1]
                technician = parts[runIdentifier_index + 2]
                batch = parts[runIdentifier_index + 3]
                run = parts[runIdentifier_index + 4]
            else:
                barcode = None
                technician = parts[runIdentifier_index + 2]
                batch = parts[runIdentifier_index + 3]
                run = parts[runIdentifier_index + 4]

        # Extract metadata from directory structure
        dir_parts = storage_location.split('/')

        # Library strategy
        if "_wgs" in dir_parts[4]:
            library_strategy = "WGS"
        else:
            library_strategy = dir_parts[5].split('_')[1]

        # Platform
        platform = "ILLUMINA" if "illumina" in dir_parts[4] else "Unknown"

        # Instrument model
        instrument_model = "Illumina_HiSeqx10" if "hiseqx10" in dir_parts[4] else "Unknown"

        # Update the NCIG_ids table in the SQLite database
        cursor.execute("""
            UPDATE NCIG_ids 
            SET Flowcell_ID = ?, 
                Lane = ?, 
                DateSequenced = ?, 
                Sample_ID = ?, 
                Species = ?, 
                Barcode = ?, 
                RunIdentifier = ?, 
                DateCollected = ?, 
                Technician = ?, 
                Batch = ?, 
                Run = ?, 
                LibraryStrategy = ?, 
                Platform = ?, 
                InstrumentModel = ?
            WHERE filename1 = ?
        """, (flowcellid, lane, date_sequenced, sampleID, species, barcode, runIdentifier, date_collected, technician, batch, run, library_strategy, platform, instrument_model, storage_location))

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
