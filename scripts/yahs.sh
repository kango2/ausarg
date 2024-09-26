#!/bin/bash
#PBS -N yahs
#PBS -P xl04
#PBS -q normal
#PBS -l walltime=01:00:00
#PBS -l mem=16GB
#PBS -l ncpus=8
#PBS -l storage=gdata/xl04+gdata/if89
#PBS -l wd
#PBS -M kirat.alreja@anu.edu.au
#PBS -j oe

module load yahs

#yahs ${ref} ${bam} -e GATC,GANTC,CTNAG,TTAA -r 10000,20000,50000,100000,200000,500000,1000000,1500000 -o ${sample}
yahs ${ref} ${bam} -e GATC,GANTC,CTNAG,TTAA -r 10000,20000,50000,100000,200000,500000,1000000,1500000 -o ${sample}