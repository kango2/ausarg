#!/bin/bash
#PBS -N juicer
#PBS -P xl04
#PBS -q normalsr
#PBS -l walltime=48:00:00
#PBS -l mem=512GB
#PBS -l ncpus=104
#PBS -l storage=gdata/xl04+gdata/if89
#PBS -l wd
#PBS -j oe
#PBS -M kirat.alreja@anu.edu.au

module load seqkit bwa python3 parallel
source /g/data/xl04/ka6418/miniconda/etc/profile.d/conda.sh
conda activate autohic
mkdir -p ${WORK_DIR}/fastq
ln -s ${FASTQ_DIR}/* ${WORK_DIR}/fastq
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
/g/data/xl04/ka6418/HiC_Assembly_Experiement/juicer/scripts/juicer.sh -D $path_juicer -d ${WORK_DIR} -y ${WORK_DIR}/base_Arima.txt -t ${PBS_NCPUS}  -g base -s Arima -z base.fa -p base.sizes >juicer.log.o 2>juicer.log.e
wait
awk -f $path_3d/utils/generate-assembly-file-from-fasta.awk base.fa >base.assembly 2>generate.log.e
wait
$path_3d/visualize/run-assembly-visualizer.sh base.assembly aligned/merged_nodups.txt >visualizer.log.o 2>visualizer.log.e
wait
python $path_greenhill/utils/fasta_to_juicebox_assembly.py base.fa >base.ctg_info.assembly

#$path_3d/run-asm-pipeline.sh --build-gapped-map -r 2 ${REF} aligned/merged_nodups.txt > 3d_dna.log 