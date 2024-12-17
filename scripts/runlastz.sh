#!/bin/bash
#PBS -q normalsr
#PBS -l ncpus=104,mem=500GB,walltime=48:00:00,jobfs=400GB,storage=gdata/te53+gdata/if89+gdata/xl04
#PBS -j oe

module load lastz/1.04.15 samtools/1.19.2

# Function to check if a file is gzipped
is_gzipped() {
    local file=$1
    testgz=$(file "$file" | grep -c 'gzip compressed data')
    if [ $testgz -eq 1 ]
    then
        return 1
    else
        return 0
    fi
}

# Function to unzip a gzipped file
unzip_gzipped() {
    local file=$1
    local outbase=$2
    pigz --keep --stdout --processes ${PBS_NCPUS} --decompress ${file} > ${PBS_JOBFS}/${outbase}
}

cd ${PBS_JOBFS}

if is_gzipped ${queryfa}
then
    unzip_gzipped ${queryfa} query.fa
else
    cp ${queryfa} query.fa
fi

if is_gzipped ${targetfa}
then
    unzip_gzipped ${targetfa} target.fa
else
    cp ${targetfa} target.fa
fi


# Index the query genome
samtools faidx query.fa

cut -f1 ${inputfasta}.fai | xargs -P ${PBS_NCPUS} -I{} bash -c 'samtools faidx query.fa {} > {}.fa'
cut -f1 ${inputfasta}.fai | xargs -P ${PBS_NCPUS} -I{} bash -c 'trf-mod {}.fa > {}.trf'

# gzip output files
