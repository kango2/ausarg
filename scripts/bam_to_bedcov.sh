#!/bin/bash
#PBS -N bedcov
#PBS -P xl04
#PBS -q normal
#PBS -l storage=gdata/xl04+gdata/if89+gdata/te53
#PBS -l walltime=2:00:00
#PBS -l mem=192GB
#PBS -l ncpus=48
#PBS -l wd
#PBS -j oe 

#usage qsub -v bam=,window=,outdir= 

module load samtools parallel
set -ex

export WINDOW=${window}
bambase="$(basename ${bam} .bam)"

samtools index -b -@ ${PBS_NCPUS} ${bam}
# Create BED file
samtools view -H ${bam} | perl -lne 'if ($_=~/SN:(\S+)\tLN:(\d+)/){ $c=$1;$l=$2; for ($i=0;$i<$l;$i+=$ENV{"WINDOW"}) { print "$c\t$i\t". ((($i+$ENV{"WINDOW"}) > $l) ? $l : ($i+$ENV{"WINDOW"}))  }} ' > "${outdir}/${bambase}.${window}.bed"

# Function to calculate depth for each BED region
calculate_depth() {
    region=$1
    samtools bedcov <(echo "$region") ${bam} | awk -v window=${WINDOW} '{print $1, $2, $3, $4/window}'
}

export -f calculate_depth

# Run depth calculation in parallel
parallel --will-cite -a "${outdir}/${bambase}.${window}.bed" -j ${PBS_NCPUS} calculate_depth > "${outdir}/${bambase}.${window}.depth.tmp.bed"
sort -k1,1 -k2,2n "${outdir}/${bambase}.${window}.depth.tmp.bed" > "${outdir}/${bambase}.${window}.depth.bed"
rm "${outdir}/${bambase}.${window}.bed"
rm "${outdir}/${bambase}.${window}.depth.tmp.bed"
touch ${outdir}/${bambase}.${window}.depth.bed.done
