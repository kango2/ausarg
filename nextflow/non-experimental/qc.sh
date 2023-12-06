#!/bin/bash
#PBS -N qc_flow
#PBS -P xl04
#PBS -q copyq
#PBS -l walltime=10:00:00
#PBS -l mem=4GB
#PBS -l ncpus=1
#PBS -l storage=gdata/xl04+gdata/if89
#PBS -l wd
#PBS -M kirat.alreja@anu.edu.au

module load nextflow
nextflow /g/data/xl04/ka6418/github/ausarg/nextflow/non-experimental/onlyqc.nf -config /g/data/xl04/ka6418/github/ausarg/nextflow/non-experimental/genejigsaw.config --topfolder /g/data/xl04/bpadata/Bassiana_duperreyi/raw/evaluation/longread_qc 