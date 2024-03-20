#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <PAF file> <UserInput1> <UserInput2> <outfile>"
    exit 1
fi

PAF_FILE=$1
USER_INPUT1=$2
USER_INPUT2=$3
OUTPUT_FILE=$4 # Specify the name of the output file

# Write the header to the output file
echo -e "Genome\tHitGenome\tSeqName\tStart\tEnd\tStrand\tHit\tHitStart\tHitEnd" > "$OUTPUT_FILE"

# Process the PAF file and append the data to the output file
while IFS=$'\t' read -r col1 col2 col3 col4 col5 col6 col7 col8 col9 col10 col11 col12; do
    echo -e "$USER_INPUT1\t$USER_INPUT2\t$col1\t$col3\t$col4\t$col5\t$col6\t$col8\t$col9" >> "$OUTPUT_FILE"
done < "$PAF_FILE"
