#!/bin/bash -ue
REF=/g/data/xl04/ka6418/github/ausarg/nextflow/outtest/yahs/*_YAHS.fasta
R1=/g/data/xl04/ka6418/nextflow_testing/testdata/fastq/HIC_R1.fastq.gz
R2=/g/data/xl04/ka6418/nextflow_testing/testdata/fastq/HIC_R2.fastq.gz
output=/g/data/xl04/ka6418/github/ausarg/nextflow/outtest/juicer


module load seqkit bwa python3 parallel
source /g/data/xl04/ka6418/miniconda/etc/profile.d/conda.sh
conda activate autohic

label=(basename ${REF} .fasta)
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
echo "Starting Juicer"
bash /g/data/xl04/ka6418/HiC_Assembly_Experiement/juicer/scripts/juicer.sh -D $path_juicer -d ${WORK_DIR} -y ${WORK_DIR}/base_Arima.txt -t ${PBS_NCPUS}  -g base -s Arima -z base.fa -p base.sizes >juicer.log.o 2>juicer.log.e
wait
echo "Juicer finished"
echo "Running 3D-DNA"
awk -f $path_3d/utils/generate-assembly-file-from-fasta.awk base.fa >base.assembly 2>generate.log.e
echo "finisheed .assembly"
wait
bash $path_3d/visualize/run-assembly-visualizer.sh base.assembly aligned/merged_nodups.txt || true
echo "finisheed .hic"
wait
python $path_greenhill/utils/fasta_to_juicebox_assembly.py base.fa >base.ctg_info.assembly
echo "finisheed final"
echo "moving things now"
mv base.hic ${label}_HifiASM_YAHS.hic
mv base.ctg_info.assembly ${label}_HifiASM_YAHS.ctg_info.assembly
mv ${label}_HifiASM_YAHS.hic ${output}
mv ${label}_HifiASM_YAHS.ctg_info.assembly ${output}
mv aligned/merged_nodups.txt ${output}/${label}_HifiASM_YAHS_JuicerAlignment.txt
echo "finished move"
