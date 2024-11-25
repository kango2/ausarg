#!/bin/bash
#PBS -N juicer
#PBS -P xl04
#PBS -q normal
#PBS -l walltime=48:00:00
#PBS -l mem=192GB
#PBS -l ncpus=48
#PBS -l storage=gdata/xl04+gdata/if89
#PBS -l wd
#PBS -j oe
#PBS -M kirat.alreja@anu.edu.au


#Sourced from GreenHill

module load seqkit bwa python3 parallel
source /g/data/xl04/ka6418/miniconda/etc/profile.d/conda.sh
conda activate autohic

PBS_JOBFS=${outputdir}
label=$(basename ${REF} .fasta)
WORK_DIR=${PBS_JOBFS}
mkdir -p ${WORK_DIR}/fastq
ln -s ${R1} ${WORK_DIR}/fastq/${label}_R1.fastq.gz
ln -s ${R2} ${WORK_DIR}/fastq/${label}_R2.fastq.gz


cd ${WORK_DIR}
seqkit sort -lr ${REF} > base.fa
bwa index base.fa > bwa_index.log 
wait
seqkit fx2tab -nl base.fa > base.sizes
JUICER_DIR=/g/data/xl04/ka6418/HiC_Assembly_Experiement/juicer
python $JUICER_DIR/misc/generate_site_positions.py Arima base ${WORK_DIR}/base.fa
wait
echo "Generate site positions done"
path_juicer=/g/data/xl04/ka6418/HiC_Assembly_Experiement/juicer
path_3d=/g/data/xl04/ka6418/HiC_Experiment/tools/3d-dna
path_greenhill=/g/data/xl04/ka6418/greenhill/GreenHill
bash /g/data/xl04/ka6418/HiC_Assembly_Experiement/juicer/scripts/juicer.sh -D $path_juicer -d ${WORK_DIR} -y ${WORK_DIR}/base_Arima.txt -t ${PBS_NCPUS}  -g base -s Arima -z base.fa -p base.sizes >juicer.log.o 2>juicer.log.e
wait
awk -f $path_3d/utils/generate-assembly-file-from-fasta.awk base.fa >base.assembly 2>generate.log.e
wait
$path_3d/visualize/run-assembly-visualizer.sh base.assembly aligned/merged_nodups.txt >visualizer.log.o 2>visualizer.log.e
wait
python $path_greenhill/utils/fasta_to_juicebox_assembly.py base.fa >base.ctg_info.assembly

#mv base.hic ${label}_Juicer.hic
#mv base.ctg_info.assembly ${label}_Juicer.ctg_info.assembly
#mv base.hic ${output}
#mv base.ctg_info.assembly ${output}
