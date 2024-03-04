#!/bin/bash
#PBS -N fast5toslow5
#PBS -P xl04
#PBS -q hugemem
#PBS -l walltime=48:00:00
#PBS -l mem=192GB
#PBS -l ncpus=48
#PBS -l jobfs=1370GB
#PBS -l storage=gdata/xl04+gdata/if89
#PBS -l wd
#PBS -M kirat.alreja@anu.edu.au

set -x

#qsub -v tarfile=,sample=,outdir= /g/data/xl04/ka6418/temp/fast5tartoslow5/sample.sh

# Change to the PBS job filesystem directory
cd ${PBS_JOBFS}
echo "Changed to directory: ${PBS_JOBFS}"

# Load slow5tools module
module load slow5tools
echo "Loaded slow5tools module"

# Unpack the tar file
echo "Unpacking tar file: ${tarfile}"
tar -xvf ${tarfile}

# Extract folder name from tar file
foldername="$(basename "${tarfile}" .tar)"
echo "Extracted folder name: ${foldername}"

# Change to the directory containing fast5 files
cd ${PBS_JOBFS}/${foldername}
echo "Changed to directory: ${PBS_JOBFS}/${foldername}"

# Decompress the fast5 files
echo "Decompressing fast5 files"
gunzip *.fast5.gz

# Create a directory for the sample
echo "Creating directory for sample: ${PBS_JOBFS}/${sample}"
mkdir -p ${PBS_JOBFS}/${sample}

# Set the blow5 directory
blow5_dir="${PBS_JOBFS}/${sample}"
echo "Set blow5 directory: ${blow5_dir}"

# Convert fast5 to blow5 format
echo "Converting fast5 to blow5 format"
slow5tools f2s ${PBS_JOBFS}/${foldername} -d ${blow5_dir}

# Merge blow5 files
echo "Merging blow5 files into ${outdir}/${sample}.blow5"
slow5tools merge -t ${PBS_NCPUS} -o "${outdir}/${sample}.blow5" ${blow5_dir}



