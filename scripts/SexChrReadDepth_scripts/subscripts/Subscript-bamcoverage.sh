#!/bin/bash
#PBS -N bamcoverage
#PBS -l ncpus=12,mem=120gb,walltime=12:00:00,storage=gdata/if89+gdata/xl04
#PBS -j oe
#PBS -M z5205618@ad.unsw.edu.au
#PBS -m ae

export input_path="/path/to/directory/containing/all/input.bam"
export binsize="20000" #20kb is good for genomes with largest chr at ~300Mbp
export workingdir="/path/to/workingdir"
export processScript="/path/to/process_bedgraph.py" #should be available in the same directory as this script on github

module use /g/data/if89/apps/modulefiles
module load parallel/20191022 pythonlib/3.9.2 python3-as-python
mkdir -p ${workingdir}/bamcoverage
mkdir -p ${workingdir}/bamcoverage/raw
mkdir -p ${workingdir}/bamcoverage/fixed

cd ${input_path}
for i in $(ls *.bam); do echo ${i/.bam/}; done | \
parallel --jobs ${PBS_NCPUS} bamCoverage -bs ${binsize} -p ${PBS_NCPUS} -of bedgraph -b {}.bam -o ${workingdir}/bamcoverage/raw/{}.bedgraph

cd ${workingdir}/bamcoverage/raw
for i in $(ls *.bedgraph); do python ${processScript} ${i} ${binsize} > ${workingdir}/bamcoverage/fixed/${i/.bedgraph/}_fixed.bedgraph; done
