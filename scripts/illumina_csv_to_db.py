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
    cursor.execute("SELECT File_path FROM Illumina_Metrics")
    existing_file_paths = set([row[0] for row in cursor.fetchall()])

    for csv_file in csv_files:
        csv_path = os.path.join(csv_folder, csv_file)

        # Read the CSV file into a DataFrame
        df = pd.read_csv(csv_path)

        # Filter rows where File_path is already in the database
        df = df[~df['File_path'].isin(existing_file_paths)]

        if not df.empty:
            # Append the non-duplicate data from the DataFrame into the SQLite table
            df.to_sql('Illumina_Metrics', conn, if_exists='append', index=False)
            
            # Update the status column to "completed" for rows where the status is "Pending"
            cursor.execute("UPDATE Illumina_Metrics SET status = 'Completed' WHERE status = 'Running'")
            conn.commit()

        # Update the set of existing File_path values for subsequent CSV checks
        existing_file_paths.update(df['File_path'].tolist())

    conn.close()
    print(f"Imported data from {len(csv_files)} CSV files to the database, skipping duplicates based on File_path.")

# Note: This function definition is meant for use outside this notebook environment, especially with argparse in use.
# If you'd like to test the function within this notebook, we'd have to call the function directly with parameters.

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Import CSV data into a SQLite database.")
    
    # Add command line arguments with flags
    parser.add_argument("-db", required=True, type=str, help="Path to the SQLite database.")
    parser.add_argument("-output", required=True, type=str, help="Path to the folder containing CSV files.")
    
    # Parse the arguments
    args = parser.parse_args()
    
    # Call the function using the parsed arguments
    import_csv_to_sqlite_v3(args.db, args.output)

