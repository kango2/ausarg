#!/bin/bash -ue
module load samtools parallel picard bwa


REF=/g/data/xl04/ka6418/nextflow_testing/testdata/BASDU_HifiASM.fasta
R1=/g/data/xl04/ka6418/nextflow_testing/testdata/fastq/HIC_R1.fastq.gz
R2=/g/data/xl04/ka6418/nextflow_testing/testdata/fastq/HIC_R2.fastq.gz
output=/g/data/xl04/ka6418/nextflow_testing/pipelinetest/scaffolding/arima

export JAVA_HOME=/apps/java/jdk-17.0.2/bin/java
label=$(basename ${REF} .fasta)
WORK_DIR=${PBS_JOBFS}

cp $REF $WORK_DIR
cd $WORK_DIR


mkdir fastq
ln -s ${R1} fastq/${label}_R1.fastq.gz
ln -s ${R2} fastq/${label}_R2.fastq.gz

IN_DIR=fastq


GENOME_BASENAME=$(basename "$REF")
REF=$WORK_DIR/$GENOME_BASENAME
SRA=$label
LABEL=$label
PREFIX=$GENOME_BASENAME
FAIDX="${GENOME_BASENAME}.fai"
RAW_DIR="$WORK_DIR/raw"
FILT_DIR="$WORK_DIR/filt"
TMP_DIR="$WORK_DIR/temp"
PAIR_DIR="$WORK_DIR/pair"
REP_DIR="$WORK_DIR/dedup"
MERGE_DIR="$WORK_DIR/merge"
REP_LABEL="${LABEL}_rep1"
FILTER="/g/data/xl04/ka6418/arima/mapping_pipeline/filter_five_end.pl"
COMBINER="/g/data/xl04/ka6418/arima/mapping_pipeline/two_read_bam_combiner.pl"
STATS="/g/data/xl04/ka6418/arima/mapping_pipeline/get_stats.pl"

MAPQ_FILTER=10
CPU=${PBS_NCPUS}

echo "### Step 0: Check output directories' existence & create them as needed"
[ -d $RAW_DIR ] || mkdir -p $RAW_DIR
[ -d $FILT_DIR ] || mkdir -p $FILT_DIR
[ -d $TMP_DIR ] || mkdir -p $TMP_DIR
[ -d $PAIR_DIR ] || mkdir -p $PAIR_DIR
[ -d $REP_DIR ] || mkdir -p $REP_DIR
[ -d $MERGE_DIR ] || mkdir -p $MERGE_DIR


if [[ -f "$REF.amb" && -f "$REF.ann" && -f "$REF.bwt" && -f "$REF.pac" && -f "$REF.sa" ]]; then
echo "BWA index exists"
else
    echo "BWA index does not exist, creating index..."
    bwa index -a bwtsw -p $PREFIX $REF
fi

if [[ -f "$REF.fai" ]]; then
    echo "Samtools index exists"
else
    echo "Samtools index does not exist, creating index..."
    samtools faidx $REF
fi

HALF_CPU=$(($PBS_NCPUS / 2))
echo $REF 
echo $GENOME_BASENAME

#Step 1: FASTQ to BAM
echo -e "${SRA}_R1.fastq.gz\n${SRA}_R2.fastq.gz" | parallel -j 2 --eta     "bwa mem -t $HALF_CPU $REF $IN_DIR/{} | samtools view -@ $HALF_CPU -Sb - > $RAW_DIR/{/.}.bam"

# Ensure the bwa mem processes have completed before moving on
wait

# Step 2: Filter 5' end
echo -e "${SRA}_R1.fastq.bam\n${SRA}_R2.fastq.bam" | parallel -j 2 --eta     "samtools view -h $RAW_DIR/{} | perl $FILTER | samtools view -Sb - > $FILT_DIR/{}"

#echo "### Step 3A: Pair reads & mapping quality filter"
perl $COMBINER $FILT_DIR/${SRA}_R1.fastq.bam $FILT_DIR/${SRA}_R2.fastq.bam samtools $MAPQ_FILTER | samtools view -bS -t $FAIDX - | samtools sort -@ $CPU -o $TMP_DIR/$SRA.bam -

#echo "### Step 3.B: Add read group"
java -jar /g/data/if89/apps/picard/2.27.4/picard.jar AddOrReplaceReadGroups INPUT=$TMP_DIR/$SRA.bam OUTPUT=$PAIR_DIR/$SRA.bam ID=$SRA LB=$SRA SM=$LABEL PL=ILLUMINA PU=none

#echo "### Step 4: Mark duplicates"
java -jar /g/data/if89/apps/picard/2.27.4/picard.jar MarkDuplicates INPUT=$PAIR_DIR/$SRA.bam OUTPUT=$REP_DIR/$REP_LABEL.bam METRICS_FILE=$REP_DIR/metrics.$REP_LABEL.txt TMP_DIR=$TMP_DIR ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=TRUE

samtools index $REP_DIR/$REP_LABEL.bam

cd dedup 
cp ${label}_rep1.bam ${label}_ArimaHiC.bam 
cp ${label}_rep1.bam.bai ${label}_ArimaHiC.bam.bai
cp ${label}_ArimaHiC.bam ${output}
cp ${label}_ArimaHiC.bam.bai ${output}
