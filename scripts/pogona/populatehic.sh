#!/bin/bash

# Define the database and CSV file
DATABASE="/g/data/xl04/bpadata/ausarg.db"
HIC_CSV="/g/data/xl04/bpadata/Pogona_vitticeps/raw/db/hic.csv"

# Function to insert data from CSV to the 'hic' table
insert_hic_data() {
    local csv_file=$1
    while IFS=, read -r sample lane file; do
        # Skip the header row
        if [ "$sample" != "sample" ]; then
            sqlite3 $DATABASE <<EOF
INSERT INTO hic (sample, lane, file) 
VALUES ('$sample', '$lane', '$file');
EOF
        fi
    done < "$csv_file"
}

# Insert data into 'hic' table
insert_hic_data "$HIC_CSV"

echo "Data has been inserted into 'hic' table from $HIC_CSV."
