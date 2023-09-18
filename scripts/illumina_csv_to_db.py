import argparse
import pandas as pd
import sqlite3
import os

def import_csv_to_sqlite_v3(db_path, csv_folder):
    # Connect to the SQLite database
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    # Detect all CSV files in the specified folder
    csv_files = [f for f in os.listdir(csv_folder) if f.endswith('.csv')]

    # Get the existing File_path values from the SQLite table for duplication check
    cursor.execute("SELECT File_path, status FROM Illumina_Metrics")
    existing_data = cursor.fetchall()

    for csv_file in csv_files:
        csv_path = os.path.join(csv_folder, csv_file)

        # Read the CSV file into a DataFrame
        df = pd.read_csv(csv_path)

        # Loop through each row in the DataFrame
        for index, row in df.iterrows():
            file_path = row['File_path']
            status = row['status']

            # Check if the file path exists in the database with "Running" status
            if (file_path, status) in existing_data and status == 'Running':
                # Replace the existing row with data from the CSV file and update status to "Completed"
                cursor.execute("DELETE FROM Illumina_Metrics WHERE File_path = ? AND status = ?", (file_path, 'Running'))
                conn.commit()
                df.loc[index, 'status'] = 'Completed'

        # Append the data from the DataFrame into the SQLite table
        df.to_sql('Illumina_Metrics', conn, if_exists='append', index=False)

    conn.close()
    print(f"Imported data from {len(csv_files)} CSV files to the database, updating 'Running' rows to 'Completed'.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Import CSV data into a SQLite database.")
    
    # Add command line arguments with flags
    parser.add_argument("-db", required=True, type=str, help="Path to the SQLite database.")
    parser.add_argument("-output", required=True, type=str, help="Path to the folder containing CSV files.")
    
    # Parse the arguments
    args = parser.parse_args()
    
    # Call the function using the parsed arguments
    import_csv_to_sqlite_v3(args.db, args.output)
