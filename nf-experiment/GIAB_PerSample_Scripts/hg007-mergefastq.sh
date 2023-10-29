#!/bin/bash
#PBS -P te53
#PBS -l ncpus=16
#PBS -l mem=32GB
#PBS -l walltime=24:00:00
#PBS -l storage=scratch/te53+gdata/te53+gdata/if89
#PBS -j oe
#PBS -o merge.log

# Define paths and directories
base_dir="/g/data/te53/referencedata/openaccess/giab/HG007_NA24695-hu38168_mother/NA24695_Mother_HiSeq100x/NA24695_Mother_HiSeq100x_fastqs"
output_dir="/g/data/te53/kh3349/Genome_Reference_Materials/merged_HG007-Mother"

# Define the list of flow cells for each group
flowcells=(
"141008_D00360_0059_BHB66DADXX" 
"141015_D00360_0061_BHB66HADXX"  
"141020_D00360_0063_BHB65FADXX" 
"141117_D00360_0066_BHB7AUADXX"
)

# Loop over each flow cell in the group
for flowcell in "${flowcells[@]}"; do
    # Loop over each lot within the flow cell
    for lot_dir in "${base_dir}/${flowcell}/Sample_NA24695/"; do
        if [ -d "$lot_dir" ]; then
            # Array to hold file paths for sorting
            files_r1=()
            files_r2=()
            
            for lane in "L001" "L002"; do
                # Collect R1 files
                for file in "${lot_dir}"/*"${lane}"_R1*.fastq.gz; do
                    files_r1+=("$file")
                done

                # Collect R2 files
                for file in "${lot_dir}"/*"${lane}"_R2*.fastq.gz; do
                    files_r2+=("$file")
                done
            done

            # Sort files based on the order of numbers after R1 and R2
            sorted_files_r1=($(printf "%s\n" "${files_r1[@]}" | sort -t '_' -k5 -n))
            sorted_files_r2=($(printf "%s\n" "${files_r2[@]}" | sort -t '_' -k5 -n))

            # Merge the R1 files
            cat "${sorted_files_r1[@]}" >> "${output_dir}/merged_${flowcell}_R1.fastq.gz"

            # Merge the R2 files
            cat "${sorted_files_r2[@]}" >> "${output_dir}/merged_${flowcell}_R2.fastq.gz"
        fi
    done
done

