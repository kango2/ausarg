#!/bin/bash
#PBS -l ncpus=16,mem=64GB,walltime=01:00:00,storage=gdata/xl04+gdata/if89+gdata/te53
#PBS -j oe
#PBS -o /g/data/xl04/ka6418/bassiana/publication-v2/logs/

module load subread/2.0.6

command="subread-align -n 30 --multiMapping -B 1 -T ${PBS_NCPUS} --sortReadsByCoordinates \
-o ${outbam} \
-i ${subreadidx} \
-t 0 \
-r ${r1fq}"

if [ -n "${r2fq}" ]; then
    command+=" -R ${r2fq}"
fi

eval $command
