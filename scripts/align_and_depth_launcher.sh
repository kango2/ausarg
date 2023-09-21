#!/bin/bash

# Read the input file
input_file="$1"

# Initialize variables
declare -A paths

# Usage:
# ./align_and_depth_launcher.sh <input_file.txt>
#
# Input File Format:
# >illumina
# path_to_illumina_reads
# >ont
# path_to_ont_reads
# >pacbio
# path_to_pacbio_reads
# >reference
# path_to_reference_genome
# >output
# path_for_output_directory

# Parse the input file
while IFS= read -r line; do
    if [[ "$line" == ">illumina" ]]; then read -r; paths["illumina"]=$REPLY; fi
    if [[ "$line" == ">ont" ]]; then read -r; paths["ont"]=$REPLY; fi
    if [[ "$line" == ">pacbio" ]]; then read -r; paths["pacbio"]=$REPLY; fi
    if [[ "$line" == ">reference" ]]; then read -r; paths["reference"]=$REPLY; fi
    if [[ "$line" == ">output" ]]; then read -r; paths["output"]=$REPLY; fi
done < "$input_file"

# Extract the base name of the reference
ref_base=$(basename "${paths["reference"]}" .fasta)  # Assuming .fasta extension, modify as needed

# Generate qsub commands for each technology
echo "qsub -o /g/data/xl04/ka6418/alignment -l storage=gdata/xl04+gdata/if89 -v platform=illumina,rawreads=${paths["illumina"]},reference=${paths["reference"]},output=${paths["output"]}/${ref_base}_illumina alignment_script.pbs" >> commands.sh
echo "qsub -o /g/data/xl04/ka6418/alignment -l storage=gdata/xl04+gdata/if89 -v platform=ont,rawreads=${paths["ont"]},reference=${paths["reference"]},output=${paths["output"]}/${ref_base}_ont alignment_script.pbs" >> commands.sh
echo "qsub -o /g/data/xl04/ka6418/alignment -l storage=gdata/xl04+gdata/if89 -v platform=pacbio,rawreads=${paths["pacbio"]},reference=${paths["reference"]},output=${paths["output"]}/${ref_base}_pacbio alignment_script.pbs" >> commands.sh
