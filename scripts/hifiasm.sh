#!/bin/bash
#PBS -N hifiasm
#PBS -P xl04
#PBS -q hugemem
#PBS -l walltime=48:00:00
#PBS -l mem=1470GB
#PBS -l ncpus=48
#PBS -l jobfs=400GB
#PBS -l storage=gdata/xl04+gdata/if89
#PBS -l wd
#PBS -M kirat.alreja@anu.edu.au
#PBS -j oe 

module load hifiasm/0.19.8

#Run HifiASM
hifiasm -o "${outputpath}/${sample}" -t ${PBS_NCPUS} --write-paf --write-ec --ul ${ont} --h1 ${r1} --h2 ${r2} ${hifi}

#Convert gfa to fasta
awk '/^S/{print ">"$2;print $3}' "${outputpath}/*.hic.p_ctg.gfa" > "${outputpath}/${sample}.fasta"


