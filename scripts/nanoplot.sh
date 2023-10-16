#!/bin/bash
#PBS -N long_read_qc
#PBS -P xl04
#PBS -q normal
#PBS -l storage=gdata/xl04+gdata/if89
#PBS -l walltime=48:00:00
#PBS -l mem=64GB
#PBS -l ncpus=8
#PBS -l wd

source /g/data/xl04/ka6418/miniconda/etc/profile.d/conda.sh
conda activate nanoplot_env

base=$(basename "$fastq" .fastq)

NanoPlot -t ${PBS_NCPUS} -o ${outputdir} --tsv_stats --raw --info_in_report --N50 --dpi 300 --prefix ${base} --format pdf --fastq ${fastq}

