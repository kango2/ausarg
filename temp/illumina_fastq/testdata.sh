#!/bin/bash
#PBS -N minimapalignment
#PBS -P xl04
#PBS -q normal
#PBS -l walltime=0:10:00
#PBS -l mem=16GB
#PBS -l ncpus=8
#PBS -l storage=gdata/xl04+gdata/if89
#PBS -l wd
#PBS -M kirat.alreja@anu.edu.au


zcat /g/data/xl04/bpadata/Tiliqua_rugosa/raw/illumina/dnaseq/350746_AusARG_UNSW_HTYH7DRXX_TCATAGATTG-CACCTTAATC_S3_L001_R1_001.fastq.gz | awk 'NR<=400' | gzip > /g/data/xl04/ka6418/ausarg/temp/alternate_test.fastq.gz