#!/bin/bash
#PBS -N hicscaffolding
#PBS -P xl04
#PBS -q normal
#PBS -l storage=gdata/xl04+gdata/if89+gdata/te53
#PBS -l walltime=48:00:00
#PBS -l mem=4GB
#PBS -l ncpus=1
#PBS -l wd
#PBS -j oe

module load nextflow


nextflow run /g/data/xl04/ka6418/github/ausarg/workflows/hic/hic.nf -config /g/data/xl04/ka6418/github/ausarg/workflows/hic/hic.config -profile NCI --inputcsv /g/data/xl04/ka6418/temp/hicwf/input/sample.csv -resume
