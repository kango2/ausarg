#!/bin/bash
#PBS -N bedcov
#PBS -P xl04
#PBS -q normal
#PBS -l storage=gdata/xl04+gdata/if89
#PBS -l walltime=48:00:00
#PBS -l mem=4GB
#PBS -l ncpus=1
#PBS -l wd



module load samtools 

export WINDOW=1000
bambase="$(basename ${bam} .bam)"

samtools view -H ${bam} | perl -lne 'if ($_=~/SN:(\S+)\tLN:(\d+)/){ $c=$1;$l=$2; for ($i=0;$i<$l;$i+=$ENV{"WINDOW"}) { print "$c\t$i\t". ((($i+$ENV{"WINDOW"}) > $l) ? $l : ($i+$ENV{"WINDOW"}))  }} ' > "${outdir}/${bambase}.bed"
samtools bedcov "${bambase}".bed ${bam} | awk '{print $1, $2, $3, $4/1000}' > "${outdir}/${bambase}.depth.bed"
