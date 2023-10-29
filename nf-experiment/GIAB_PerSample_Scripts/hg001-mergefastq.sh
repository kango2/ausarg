#!/bin/bash
#PBS -P te53
#PBS -l ncpus=16
#PBS -l mem=32GB
#PBS -l walltime=48:00:00
#PBS -l storage=scratch/te53+gdata/te53+gdata/if89
#PBS -j oe
#PBS -o merge.log

# Define paths and directories
base_dir="/g/data/te53/referencedata/openaccess/giab/NA12878/NIST_NA12878_HG001_HiSeq_300x"
output_dir="/g/data/te53/kh3349/Genome_Reference_Materials/merged_HG001"

# Define the list of flow cells for each group
flowcells_g1=(
  "131219_D00360_005_BH814YADXX"
  "131219_D00360_006_AH81VLADXX"
  "131223_D00360_007_BH88WKADXX"
  "131223_D00360_008_AH88U0ADXX"
  "140313_D00360_0014_AH8GGVADXX"
  "140313_D00360_0015_BH9258ADXX"
  "140407_D00360_0016_AH948VADXX"
  "140407_D00360_0017_BH947YADXX"
)

flowcells_g2=(
  "140115_D00360_0009_AH8962ADXX" 
  "140115_D00360_0010_BH894YADXX" 
  "140127_D00360_0011_AHGV6ADXX" 
  "140127_D00360_0012_BH8GVUADXX" 
  "140207_D00360_0013_AH8G92ADXX"
)    

# Loop over each group of samples
for group in "g1" "g2"; do
  if [ "$group" == "g1" ]; then
    flowcells=("${flowcells_g1[@]}")
  elif [ "$group" == "g2" ]; then
    flowcells=("${flowcells_g2[@]}")
  fi
  
  # Loop over each flow cell in the group
  for flowcell in "${flowcells[@]}"; do
    # Loop over each lot within the flow cell
    for lot_dir in "${base_dir}/${flowcell}/Project_RM8398/Sample_"*"/"; do
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
        cat "${sorted_files_r1[@]}" >> "${output_dir}/merged_${flowcell}_checked_R1.fastq.gz"

        # Merge the sorted R2 files
        cat "${sorted_files_r2[@]}" >> "${output_dir}/merged_${flowcell}_checked_R2.fastq.gz"
      fi
    done
  done
done
