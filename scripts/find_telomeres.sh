#!/bin/bash
#PBS -N Telomere
#PBS -P xl04
#PBS -q express
#PBS -l walltime=0:30:00
#PBS -l mem=16GB
#PBS -l ncpus=48
#PBS -l storage=gdata/xl04+gdata/if89
#PBS -l wd
#PBS -l jobfs=400GB
#PBS -M kirat.alreja@anu.edu.au

#module load kentutils/0.0 TRF biopython parallel

inputfile=/g/data/xl04/ka6418/tiliqua/verkko_assembly/asm/assembly.unassigned.fasta
outputdir=/g/data/xl04/ka6418/temp
percentage_match="$permatch"
number_copies="$copies"

PBS_JOBFS=/g/data/xl04/ka6418/temp/testparallel

#faSplit sequence $inputfile 1000 ${PBS_JOBFS}/chunk

cd ${PBS_JOBFS}
#filelist=$(ls ${PBS_JOBFS}/chunk*)
#printf "%s\n" "${filelist[@]}" | parallel -I{} --jobs ${PBS_NCPUS} trf {} 2 5 7 80 10 50 2000 -l 10 -d -h 


#for file in ${PBS_JOBFS}/*.dat;
#do
#    python3 /g/data/xl04/ka6418/ausarg/scripts/trf2gff.py -i ${file}
#done


#cat *.gff3 > $(basename "$inputfile" .fasta).gff3

# Create the CSV header
echo "Sequence_ID,Start,End,ID,period,copies,consensus_size,perc_match,perc_indels,align_score,entropy,cons_seq,repeat_seq" > $(basename "$inputfile" .fasta).csv

# Process GFF3 data and convert to CSV using awk
awk -F'\t' 'BEGIN {OFS=","}
    {
        split($9, attributes, /[;=]/)
        print $1, $4, $5, attributes[2], attributes[4], attributes[6], attributes[8], attributes[10], attributes[12], attributes[14], attributes[16], attributes[18], attributes[20]
    }
' $(basename "$inputfile" .fasta).gff3 >> $(basename "$inputfile" .fasta).csv

awk -F',' 'NR == 1 || ($5 == 6 && $6 > 100 && $8 >= 90)' $(basename "$inputfile" .fasta).csv  > $(basename "$inputfile" .fasta).csv
awk -F',' -v val1="$number_copies" -v val2="$percentage_match" 'NR == 1 || ($5 == 6 && $6 > val1 && $8 >= val2)' "$(basename "$inputfile" .fasta).csv" > "$(basename "$inputfile" .fasta)_filtered.csv"

#awk 'NR == 1; NR > 1 && FNR > 1' ${PBS_JOBFS}/trf/*.csv > ${PBS_JOBFS}/trf/combined.csv
#rsync -a ${PBS_JOBFS}/combined.csv ${outputdir}












