#!/bin/bash

# Define the database file
DATABASE="/g/data/xl04/bpadata/ausarg.db"

# Step 1: Add the new 'chemistry' column to the 'ont' table
sqlite3 $DATABASE <<EOF
ALTER TABLE ont ADD COLUMN chemistry TEXT;
EOF

# Step 2: Update the 'chemistry' column with the value 'R9' for all existing rows
sqlite3 $DATABASE <<EOF
UPDATE ont SET chemistry = 'R9';
EOF

# Output success message
echo "Column 'chemistry' has been added to the 'ont' table, and all rows have been updated with the value 'R9'."
