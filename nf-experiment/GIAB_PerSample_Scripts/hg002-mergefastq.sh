#!/bin/bash
#PBS -P te53
#PBS -l ncpus=16
#PBS -l mem=32GB
#PBS -l walltime=48:00:00
#PBS -l storage=scratch/te53+gdata/te53+gdata/if89
#PBS -j oe
#PBS -o merge.log

# Define paths and directories
base_dir="/g/data/te53/referencedata/openaccess/giab/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/HG002_HiSeq300x_fastq"
output_dir="/g/data/te53/kh3349/Genome_Reference_Materials/merged_HG002-Son"

# Define the list of flow cells for each group
flowcells=(
"140528_D00360_0018_AH8VC6ADXX" 
"140528_D00360_0019_BH8VDAADXX" 
"140605_D00360_0020_AH9V1RADXX" 
"140605_D00360_0021_BH9V1VADXX" 
"140609_D00360_0022_AH9UJNADXX" 
"140609_D00360_0023_BH9UD5ADXX" 
"140611_D00360_0024_AH9TKDADXX" 
"140611_D00360_0025_BH9UD6ADXX" 
"140613_D00360_0026_AHA2RRADXX" 
"140613_D00360_0027_BHA2TEADXX" 
"140616_D00360_0028_AHA2RLADXX" 
"140616_D00360_0029_BHA2WPADXX"
)

# Loop over each flow cell in the group
for flowcell in "${flowcells[@]}"; do
    # Loop over each lot within the flow cell
    for lot_dir in "${base_dir}/${flowcell}/Project_RM8391_RM8392/Sample_"*"/"; do
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
