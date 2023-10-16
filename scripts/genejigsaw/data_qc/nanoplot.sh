#!/bin/bash
#PBS -N alignment
#PBS -P xl04
#PBS -q normal
#PBS -l storage=gdata/xl04+gdata/if89
#PBS -l walltime=05:00:00
#PBS -l mem=192GB
#PBS -l ncpus=48
#PBS -l wd

source /g/data/xl04/ka6418/miniconda/etc/profile.d/conda.sh
conda activate nanoplot_env

NanoPlot -t 48 --huge -o /g/data/xl04/ka6418/genejigsaw_testing --tsv_stats --fastq /g/data/xl04/bpadata/Lampropholis_delicata/raw/pacbio/350750_AusARG_AGRF_DA064169.ccs.fq.gz

