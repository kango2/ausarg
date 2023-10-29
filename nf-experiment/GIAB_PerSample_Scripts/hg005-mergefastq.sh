#!/bin/bash
#PBS -P te53
#PBS -l ncpus=16
#PBS -l mem=32GB
#PBS -l walltime=24:00:00
#PBS -l storage=scratch/te53+gdata/te53+gdata/if89
#PBS -j oe
#PBS -o merge.log

# Define paths and directories
base_dir="/g/data/te53/referencedata/openaccess/giab/HG005_NA24631_son/HG005_NA24631_son_HiSeq_300x/basespace_250bps_fastqs"
output_dir="/g/data/te53/kh3349/Genome_Reference_Materials/merged_HG005-Son"

# Make sure the output directory exists
mkdir -p "$output_dir"

# Define the list of flow cells for each group
flowcells=(
"150420_HG005_Homogeneity_01-22889870"
"150424_HG005_Homogeneity_02_FCA-22108087"
"150424_HG005_Homogeneity_02_FCB-22108088"
"150430_HG005_Homogeneity_03_FCA-22320298"
"150430_HG005_Homogeneity_03_FCB-22320298"
"150506_HG005_Homogeneity_04_FCA-22365346"
"150506_HG005_Homogeneity_04_FCB-22365345"
)

# Loop over each flow cell in the group
for flowcell in "${flowcells[@]}"; do
    echo "Processing flow cell $flowcell..."
    # Loop over each lot within the flow cell
    for lot_dir in "${base_dir}/${flowcell}/5"*"/"; do
        if [ -d "$lot_dir" ]; then
            # Array to hold file paths for sorting
            files_r1=()
            files_r2=()

            for lane in "L001" "L002"; do
                # Collect R1 files
                for file in "${lot_dir}"*"${lane}"_R1*.fastq.gz; do
                    files_r1+=("$file")
                done
                # Collect R2 files
                for file in "${lot_dir}"*"${lane}"_R2*.fastq.gz; do
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
