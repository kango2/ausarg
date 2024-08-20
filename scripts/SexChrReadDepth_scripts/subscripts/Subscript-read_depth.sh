#!/bin/bash
#PBS -N JUST_subread
#PBS -l ncpus=12,mem=120gb,storage=gdata/if89+gdata/xl04,walltime=12:00:00
#PBS -j oe
#PBS -M z5205618@ad.unsw.edu.au
#PBS -m ae

export workingdir="/path/to/workingdir"
export genome="/path/to/genome.fasta"
export inputcsv="/path/to/inputs.csv" #inputs.csv should be a 3 column csv file with column1 as identifier, column2 as path to forward reads, column3 as path to reverse reads, an example is provided

module use /g/data/if89/apps/modulefiles
module load subread/2.0.6 parallel/20191022 samtools/1.18
cd ${workingdir}

#index genome
mkdir -p ${workingdir}/subread/index
subread-buildindex -o ${workingdir}/subread/index/genome ${genome}

#subread mapping to genome
mkdir -p ${workingdir}/subread/mapping

while IFS=',' read -r col1 col2 col3
do
    echo -e "subread-align -t 1 -T ${PBS_NCPUS} --sortReadsByCoordinates -i ${workingdir}/subread/index/genome -r ${col2} -R ${col3} -o ${workingdir}/subread/mapping/${col1}_subread.bam"
done < ${inputcsv} | parallel --jobs ${PBS_NCPUS} {}

# re-generate the .bai index files because bamcoverage can not read subread-generated bai for whatever reason
cd ${workingdir}/subread/mapping
rm *.bai
for i in $(ls *.bam); do echo ${i}; done | \
parallel --jobs ${PBS_NCPUS} samtools index -b -@ ${PBS_NCPUS} {} {}.bai
