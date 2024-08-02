#!/bin/bash
#PBS -N hifiasm
#PBS -P xl04
#PBS -q hugemem
#PBS -l walltime=48:00:00
#PBS -l mem=1470GB
#PBS -l ncpus=48
#PBS -l jobfs=400GB
#PBS -l storage=gdata/if89+gdata/xl04
#PBS -l wd
#PBS -j oe 

set -ex
set -o pipefail
set -u

usage() {
        echo "Usage: qsub -o {logs directory} -v outputdir={out directory}, sample={sample}, ont={: separated ont files}, hifi={: separated hifi files}" >&2
        echo
        exit 1
}

module load hifiasm

mkdir -p ${outputdir}

ont="${ont//:/,}"
hifi="${hifi//:/ }"

hifiasm --telo-m CCCTAA -o "${outputdir}/${sample}" -t ${PBS_NCPUS} --ul ${ont} --h1 ${h1} --h2 ${h2} ${hifi}



