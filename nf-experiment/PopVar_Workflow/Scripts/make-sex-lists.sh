#!/bin/bash

# Input list file
FILE=$1

# Empty the lists if they already exist
> female.list
> male.list

# Populate the lists based on file naming
while IFS= read -r line; do
    if [[ $line =~ Female\.raw_variants\.g\.vcf\.gz$ ]]; then
        echo $line >> female.list
    elif [[ $line =~ Male\.raw_variants\.g\.vcf\.gz$ ]]; then
        echo $line >> male.list
    fi
done < "$FILE"

