import sqlite3
import pandas as pd
import argparse

def main(db_path, sra_data_path):
    # Create a new SQLite database
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    # Read the updated "SRA_data" sheet from the provided Excel file
    updated_sra_data = pd.read_excel(sra_data_path, sheet_name="SRA_data")

    # Insert the data from the "SRA_data" sheet into the Sample table
    updated_sra_data.to_sql('SRA', conn, if_exists='replace', index=False)

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
