import sqlite3
import pandas as pd
import argparse
import sqlite3

def get_illumina_filenames_from_sra(db_path):
    """Retrieve filenames from the SRA table."""
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    cursor.execute("SELECT filename FROM SRA where library_strategy== 'WGS' AND library_source=='GENOMIC' AND platform=='ILLUMINA'")
    filenames = cursor.fetchall()
    conn.close()
    
    # Flatten the filenames list and remove None entries
    flat_filenames = [f for sublist in filenames for f in sublist if f]
    
    return flat_filenames

def populate_illumina_metrics(db_path, filenames):
    """Insert filenames into the illumina_metrics table with status set to 'Pending', avoiding duplicates using INSERT OR IGNORE."""
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    
    rows_inserted = 0  # Counter for rows inserted
    
    for filename in filenames:
        cursor.execute("INSERT OR IGNORE INTO illumina_metrics (File_path, status) VALUES (?, 'Pending')", (filename,))
        
        # If a row was inserted, the cursor's rowcount attribute will be 1. Otherwise, it'll be -1.
        if cursor.rowcount == 1:
            rows_inserted += 1
            
    conn.commit()
    conn.close()
    
    print(f"{rows_inserted} rows inserted into illumina_metrics.")

def main(db_path, sra_data_path):
    # Create a new SQLite database
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    # Read the updated "SRA_data" sheet from the provided Excel file
    updated_sra_data = pd.read_excel(sra_data_path, sheet_name="SRA_data")

    # Insert the data from the "SRA_data" sheet into the Sample table
    updated_sra_data.to_sql('SRA', conn, if_exists='replace', index=False)
    
    #Initialise the illumina metrics table
    populate_illumina_metrics(db_path,get_illumina_filenames_from_sra(db_path))

    # Close the database connection
    conn.close()

    # Confirm that data insertion is complete
    print("Data insertion into the SRA table is complete.")

    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Insert data from an Excel file into an SQLite database.")
    
    # Add flags for database and excel file
    parser.add_argument("--db", required=True, help="Path to the SQLite database.")
    parser.add_argument("--excel", required=True, help="Path to the Excel file containing SRA data.")

    args = parser.parse_args()

    main(args.db, args.excel)
