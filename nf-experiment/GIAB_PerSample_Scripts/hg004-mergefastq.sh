#!/bin/bash
#PBS -P te53
#PBS -l ncpus=16
#PBS -l mem=32GB
#PBS -l walltime=48:00:00
#PBS -l storage=scratch/te53+gdata/te53+gdata/if89
#PBS -j oe
#PBS -o merge.log

# Define paths and directories
base_dir="/g/data/te53/referencedata/openaccess/giab/HG004_NA24143_mother/NIST_HiSeq_HG004_Homogeneity-14572558/HG004_HiSeq300x_fastq"
output_dir="/g/data/te53/kh3349/Genome_Reference_Materials/merged_HG004-Mother"

# Define the list of flow cells for each group
flowcells=(
"140818_D00360_0046_AHA5R5ADXX" 
"140818_D00360_0047_BHA66FADXX" 
"140821_D00360_0048_AHAC69ADXX" 
"140821_D00360_0049_BHAC63ADXX" 
"140902_D00360_0050_AHAK0GADXX" 
"140902_D00360_0051_BHAK1WADXX" 
"140905_D00360_0052_AHAJW4ADXX" 
"140905_D00360_0053_BHAJVYADXX" 
"140915_D00360_0054_AHAJR6ADXX" 
"140915_D00360_0055_BHAML8ADXX" 
"140918_D00360_0056_AHAMJJADXX" 
"140918_D00360_0057_BHAML6ADXX"
)

# Loop over each flow cell in the group
for flowcell in "${flowcells[@]}"; do
    # Loop over each lot within the flow cell
    for lot_dir in "${base_dir}/${flowcell}/Project_RM8392/Sample_"*"/"; do
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

            # Merge the sorted R1 files
            cat "${sorted_files_r1[@]}" >> "${output_dir}/merged_${flowcell}_R1.fastq.gz"

            # Merge the sorted R2 files
            cat "${sorted_files_r2[@]}" >> "${output_dir}/merged_${flowcell}_R2.fastq.gz"
        fi
    done
done
