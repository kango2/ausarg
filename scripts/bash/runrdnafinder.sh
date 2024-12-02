#!/bin/bash
#PBS -P xl04
#PBS -q normal
#PBS -l walltime=00:20:00
#PBS -l mem=32GB
#PBS -l ncpus=8
#PBS -j oe
#PBS -o /g/data/xl04/genomeprojects/Pogona_vitticeps/logs/

module load minimap2 samtools pythonlib

# Check if the task has already been completed
checkpoint_file="${outputdir}/${sampleid}.ribocop"

if [ -f "$checkpoint_file.done" ];
then
    echo "Task already completed. Exiting."
    exit 0
else
    touch "${checkpoint_file}.running"
fi

set -ex
# Define the deuterostomia rDNA file as a variable
rdnalibfa="/g/data/xl04/hrp561/rdnalib/deuterostomia.18s28s.rDNA.fasta"

# Define the sample id as a variable
# sampleid="POGVITdef.h1.corrected_YAHS"

# Define the input fasta file as a variable
# inputfasta="POGVITdef.h1.corrected_YAHS.fasta"

# Define the output directory as a variable
# outputdir="/g/data/xl04/genomeprojects/Pogona_vitticeps/tmp/rdna/"

mkdir -p $outputdir
cd $outputdir

# Run minimap2 to generate primary PAF file
minimap2 -t ${PBS_NCPUS} --secondary=no -o ${sampleid}.primary.paf ${rdnalibfa} ${inputfasta}

# Extract the most frequent 18S and 28S sequences
eighteen=$(awk '$11>=$7*0.9' ${sampleid}.primary.paf | cut -f6 | sort | uniq -c | grep _18S | sort -k1,1n | tail -n 1 | sed -r 's/^\s+//g' | cut -f2 -d ' ')
twoeight=$(awk '$11>=$7*0.9' ${sampleid}.primary.paf | cut -f6 | sort | uniq -c | grep _28S | sort -k1,1n | tail -n 1 | sed -r 's/^\s+//g' | cut -f2 -d ' ')

# Extract sequences from the fasta file
samtools faidx $rdnalibfa $eighteen > ${sampleid}.ssulsurdna.fa
samtools faidx $rdnalibfa $twoeight >> ${sampleid}.ssulsurdna.fa

# Run minimap2 to generate refined PAF file
minimap2 -t ${PBS_NCPUS} --secondary=no -o ${sampleid}.refined.paf ${sampleid}.ssulsurdna.fa ${inputfasta}

# Run the rdnabuilder.py script with the refined PAF file
python3 /g/data/te53/hrp561/ausarg/scripts/python/rdnabuilder.py -p ${sampleid}.refined.paf -f ${inputfasta} -o ./ -s ${sampleid}

# Map assembly to the rDNA morphs
#minimap2 -t ${PBS_NCPUS} -o --secondary=no -a ${sampleid}.ssulsurdna.fa ${inputfasta} 

# TODO:
# - create a contained repo for this. perhaps include these commands in rdnabuilder.py script
# - pangenome graph for rDNA morphs
# - align assembly to morphs to ensure that all morphs are captured
# - annotate rDNA regions with repeatmasker
# - create reporting stats such as 
#   - length distribution, 
#   - completeness of morphs,
#   - head-to-tail organisation,
#   - number of scaffolds with rDNA seq only and their names
#   - number of rDNA copies in the assembly
#   - sequence divergence metrics between morphs
#   - multiple sequence alignments
#   - phylogenetic trees

# Mark the task as completed
touch "${checkpoint_file}.done"
rm -f "${checkpoint_file}.running"