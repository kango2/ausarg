#!/bin/bash
#PBS -l ncpus=48,mem=190GB,walltime=48:00:00,storage=gdata/if89+gdata/xl04,jobfs=400GB
#PBS -P xl04
#PBS -N repeatmasker
#PBS -j oe

set -ex
module use /g/data/if89/apps/modulefiles
module load RepeatMasker/4.1.2-p1

mkdir -p ${PBS_JOBFS}/rmout
mkdir -p ${outputdir}
rsync -a $inputgenome ${PBS_JOBFS}/
cd ${PBS_JOBFS}
RepeatMasker -engine rmblast -parallel ${PBS_NCPUS} -frag 1000000 -lib ${rmlib} -dir ${PBS_JOBFS}/rmout -xsmall -gff ${PBS_JOBFS}/$(basename $inputgenome)
rsync -a ${PBS_JOBFS}/rmout/ ${outputdir}/