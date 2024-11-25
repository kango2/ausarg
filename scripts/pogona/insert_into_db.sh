#!/bin/bash

# Define the database file
DATABASE="/g/data/xl04/bpadata/ausarg.db"

# CSV files
ONT_CSV="/g/data/xl04/bpadata/Pogona_vitticeps/raw/db/ont.csv"
PACBIO_CSV="/g/data/xl04/bpadata/Pogona_vitticeps/raw/db/pacbio.csv"

# Function to insert data from CSV to the database table
insert_csv_to_table() {
    local csv_file=$1
    local table_name=$2

    # Use sqlite3 to insert data into the table
    while IFS=, read -r sample flowcell file path; do
        # Skip the header row
        if [ "$sample" != "sample" ]; then
            sqlite3 $DATABASE <<EOF
INSERT INTO $table_name (sample, flowcell, file, path) 
VALUES ('$sample', '$flowcell', '$file', '$path');
EOF
        fi
    done < "$csv_file"
}

# Insert data into 'ont' table
insert_csv_to_table "$ONT_CSV" "ont"

# Insert data into 'pacbio' table
insert_csv_to_table "$PACBIO_CSV" "pacbio"

# Output success message
echo "Data has been inserted into 'ont' and 'pacbio' tables."
