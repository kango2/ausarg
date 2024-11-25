#!/bin/bash

# Define the database file
DATABASE="/g/data/xl04/bpadata/ausarg.db"

# Create the 'hic' table with the columns: sample, lane, file
sqlite3 $DATABASE <<EOF
CREATE TABLE IF NOT EXISTS hic (
    sample TEXT,
    lane TEXT,
    file TEXT
);
EOF

echo "Table 'hic' created successfully in the database."
