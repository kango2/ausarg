#!/bin/bash
#PBS -N Telomere
#PBS -q normal
#PBS -P xl04
#PBS -l storage=gdata/if89+gdata/xl04+gdata/te53
#PBS -l walltime=24:00:00
#PBS -l mem=192GB
#PBS -l ncpus=48
#PBS -l wd
#PBS -j oe
#PBS -l jobfs=100GB

#set -ex breaks the script - probably due to trf2gff conversion 

usage() {
	echo "Usage: qsub -l storage=gdata/if89+gdata/projectcode -o /path/to/stdouterr -P projectcode -v input=/path/to/fasta,output=/path/to/output/csv,permatch=90,copies=100 ./find_telomeres.sh" >&2
	echo
	exit 1
}

[ -z "${input}" ] && usage
[ -z "${output}" ] && usage
[ -z "${permatch}" ] && usage
[ -z "${copies}" ] && usage

module load kentutils/0.0 TRF/4.09.1 biopython/1.79 parallel/20191022 

#source /g/data/xl04/ka6418/miniconda/etc/profile.d/conda.sh
#conda activate autohic

inputfile="$input"
outputdir="$output"
percentage_match="$permatch"
number_copies="$copies"

faSplit sequence $inputfile 10000 ${PBS_JOBFS}/chunk

cd ${PBS_JOBFS}

num_jobs=$(($PBS_NCPUS / 4))
filelist=$(ls ${PBS_JOBFS}/chunk*)
printf "%s\n" "${filelist[@]}" | parallel -I{} --jobs ${num_jobs} trf {} 2 7 7 80 10 500 6 -l 10 -d -h 

for file in ${PBS_JOBFS}/*.dat;
do
    python3 /g/data/xl04/ka6418/github/ausarg/scripts/trf2gff.py -i ${file} -o ${file}.gff3 
done

cat *.gff3 > $(basename "$inputfile" .fasta).gff3

echo "Sequence_ID,Start,End,ID,period,copies,consensus_size,perc_match,perc_indels,align_score,entropy,cons_seq,repeat_seq" > $(basename "$inputfile" .fasta).csv

awk -F'\t' 'BEGIN {OFS=","}
    {
        split($9, attributes, /[;=]/)
        print $1, $4, $5, attributes[2], attributes[4], attributes[6], attributes[8], attributes[10], attributes[12], attributes[14], attributes[16], attributes[18], attributes[20]
    }
' $(basename "$inputfile" .fasta).gff3 >> $(basename "$inputfile" .fasta).csv

python3 /g/data/xl04/ka6418/github/ausarg/scripts/processtrftelo.py "$(basename "$inputfile" .fasta).csv" "$inputfile" "$outputdir/$(basename "$inputfile" .fasta)_Telomeres.csv" $number_copies $percentage_match












