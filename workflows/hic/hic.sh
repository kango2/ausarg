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

#Required variables 
#workdir, wfname, csv

cd ${workdir}

rm -f "${workdir}/${wfname}_report.html" "${workdir}/${wfname}_trace.txt" "${workdir}/${wfname}_timeline.html" "${workdir}/${wfname}_graph.dot"

nextflow run /g/data/xl04/ka6418/github/ausarg/workflows/hic/hic.nf \
 -work-dir ${workdir} -config /g/data/xl04/ka6418/github/ausarg/workflows/hic/hic.config \
 -profile NCI --inputcsv ${csv} -resume \
 -with-report "${workdir}/${wfname}_report.html" \
 -with-trace "${workdir}/${wfname}_trace.txt"  \
 -with-timeline "${workdir}/${wfname}_timeline.html" \
 -with-dag "${workdir}/${wfname}_graph.dot" 
 
