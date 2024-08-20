#!/bin/bash
#PBS -N blast
#PBS -l ncpus=16,walltime=6:00:00,storage=gdata/if89+gdata/xl04,mem=180GB
#PBS -j oe
#PBS -M z5205618@ad.unsw.edu.au
#PBS -m ae

export workingdir="/path/to/workingdir"
export genome="/path/to/genome.fasta"
export query="/path/to/query.fasta" #this should be the Y enriched kmers
export output="/path/to/output.tsv"


module use /g/data/if89/apps/modulefiles
module load blast/2.14.1

cd ${workingdir}
mkdir -p ${workingdir}/blast
mkdir -p ${workingdir}/blast/database

makeblastdb -in ${genome} \
-parse_seqids -title "genome" -dbtype nucl \
-out ${workingdir}/blast/database/genome

blastn -query ${query} \
-db ${workingdir}/blast/database/genome -parse_deflines -outfmt 6 -num_threads ${PBS_NCPUS} \
-max_hsps 1 -max_target_seqs 1 -out ${output}
