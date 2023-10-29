import sqlite3

database_path = "/g/data/xy86/giab/giabtest.db"

# Connect to the SQLite database
conn = sqlite3.connect(database_path)
cursor = conn.cursor()

# Fetch the required columns
query = "SELECT fastq_file1, fastq_file2, sex FROM giab_study"
cursor.execute(query)
rows = cursor.fetchall()

# Print rows in TSV format
for row in rows:
    print(f"{row[0]}\t{row[1]}\t{row[2]}")

conn.close()
