#!/bin/bash
#PBS -N qc_flow
#PBS -P xl04
#PBS -q express
#PBS -l walltime=00:30:00
#PBS -l mem=192GB
#PBS -l ncpus=48
#PBS -l storage=gdata/xl04+gdata/if89
#PBS -l wd
#PBS -M kirat.alreja@anu.edu.au


/g/data/xl04/ka6418/github/ausarg/nextflow/kmer_nf.sh -i "/g/data/xl04/ka6418/nextflow_testing/testhifi.fq.gz" -s BASDU -o "/g/data/xl04/ka6418/github/ausarg/nextflow" -l 17 -t PacBio
