#!/bin/bash

# Define the database file (you can change the name as needed)
DATABASE="/g/data/xl04/bpadata/ausarg.db"

# Remove the database if it already exists (optional)
if [ -f "$DATABASE" ]; then
    rm "$DATABASE"
fi

# Create the database and tables using sqlite3
sqlite3 $DATABASE <<EOF
-- Create 'ont' table
CREATE TABLE ont (
    sample TEXT,
    flowcell TEXT,
    file TEXT,
    path TEXT
);

-- Create 'pacbio' table
CREATE TABLE pacbio (
    sample TEXT,
    flowcell TEXT,
    file TEXT,
    path TEXT
);

-- Confirm that tables were created
.tables
EOF

# Output success message
echo "Tables 'ont' and 'pacbio' have been created in the database '$DATABASE'."
