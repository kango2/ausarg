#!/bin/bash
#PBS -P te53
#PBS -l ncpus=16
#PBS -l mem=32GB
#PBS -l walltime=48:00:00
#PBS -l storage=scratch/te53+gdata/te53+gdata/if89
#PBS -j oe
#PBS -o merge.log

# Define paths and directories
base_dir="/g/data/te53/referencedata/openaccess/giab/HG003_NA24149_father/NIST_HiSeq_HG003_Homogeneity-12389378/HG003_HiSeq300x_fastq"
output_dir="/g/data/te53/kh3349/Genome_Reference_Materials/merged_HG003-Father"

# Define the list of flow cells for each group
flowcells=(
"140627_D00360_0030_AHA0L6ADXX" 
"140701_D00360_0032_AHA0KGADXX" 
"140701_D00360_0033_BH9YY4ADXX" 
"140703_D00360_0034_AHA3FRADXX" 
"140703_D00360_0035_BHA3CGADXX" 
"140709_D00360_0036_AHA3D6ADXX" 
"140709_D00360_0037_BHA3HMADXX" 
"140711_D00360_0038_AHA3J4ADXX" 
"140711_D00360_0039_BHA3HYADXX" 
"140716_D00360_0040_AHA6C8ADXX" 
"140716_D00360_0041_BHA658ADXX" 
"140718_D00360_0042_AHA660ADXX" 
"140718_D00360_0043_BHA65JADXX" 
"140721_D00360_0044_AHA66RADXX"
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
