Channel
    .fromPath(params.inputcsv)
    .splitCsv(header: true)
    .map { row ->
        // Parse the fields
        def sample = row.sample
        def hap = row.hap
        def fasta = row.fasta
        def R1 = row.r1
        def R2 = row.r2
        def outdir = row.outputdir
        
        // Create the tuple in the desired format
        tuple("${sample}.${hap}", fasta, R1, R2, outdir)
    }
    .set {inputCh}


process arima {

    publishDir "${outputdir}", mode: 'copy'

    input: 
    tuple val(sample), val(fasta), val(R1), val(R2), val(outputdir)

    output:
    tuple val(sample), val(fasta), val(R1), val(R2), path("${sample}_ArimaHiC.bam"), val(outputdir)

    script:
    """
    module load samtools parallel picard bwa

    REF=${fasta}
    R1=${R1}
    R2=${R2}
    output=\${PWD}

    export JAVA_HOME=/apps/java/jdk-17.0.2/bin/java
    label=${sample}
    WORK_DIR=\${PBS_JOBFS}

    cp \$REF \$WORK_DIR
    cd \$WORK_DIR

    mkdir fastq
    ln -s \${R1} fastq/\${label}_R1.fastq.gz
    ln -s \${R2} fastq/\${label}_R2.fastq.gz

    IN_DIR=fastq

    GENOME_BASENAME=\$(basename "\$REF")
    REF=\$WORK_DIR/\$GENOME_BASENAME
    SRA=\$label
    LABEL=\$label
    PREFIX=\$GENOME_BASENAME
    FAIDX="\${GENOME_BASENAME}.fai"
    RAW_DIR="\$WORK_DIR/raw"
    FILT_DIR="\$WORK_DIR/filt"
    TMP_DIR="\$WORK_DIR/temp"
    PAIR_DIR="\$WORK_DIR/pair"
    REP_DIR="\$WORK_DIR/dedup"
    MERGE_DIR="\$WORK_DIR/merge"
    REP_LABEL="\${LABEL}_rep1"
    FILTER="/g/data/xl04/ka6418/arima/mapping_pipeline/filter_five_end.pl"
    COMBINER="/g/data/xl04/ka6418/arima/mapping_pipeline/two_read_bam_combiner.pl"
    STATS="/g/data/xl04/ka6418/arima/mapping_pipeline/get_stats.pl"

    MAPQ_FILTER=10
    CPU=\${PBS_NCPUS}

    echo "### Step 0: Check output directories' existence & create them as needed"
    [ -d \$RAW_DIR ] || mkdir -p \$RAW_DIR
    [ -d \$FILT_DIR ] || mkdir -p \$FILT_DIR
    [ -d \$TMP_DIR ] || mkdir -p \$TMP_DIR
    [ -d \$PAIR_DIR ] || mkdir -p \$PAIR_DIR
    [ -d \$REP_DIR ] || mkdir -p \$REP_DIR
    [ -d \$MERGE_DIR ] || mkdir -p \$MERGE_DIR


    if [[ -f "\$REF.amb" && -f "\$REF.ann" && -f "\$REF.bwt" && -f "\$REF.pac" && -f "\$REF.sa" ]]; then
    echo "BWA index exists"
    else
        echo "BWA index does not exist, creating index..."
        bwa index -a bwtsw -p \$PREFIX \$REF
    fi

    if [[ -f "\$REF.fai" ]]; then
        echo "Samtools index exists"
    else
        echo "Samtools index does not exist, creating index..."
        samtools faidx \$REF
    fi

    HALF_CPU=\$((\$PBS_NCPUS / 2))
    echo \$REF 
    echo \$GENOME_BASENAME

    #Step 1: FASTQ to BAM
    echo -e "\${SRA}_R1.fastq.gz\\n\${SRA}_R2.fastq.gz" | parallel -j 2 --eta \
    "bwa mem -t \$HALF_CPU \$REF \$IN_DIR/{} | samtools view -@ \$HALF_CPU -Sb - > \$RAW_DIR/{/.}.bam"

    # Ensure the bwa mem processes have completed before moving on
    wait

    # Step 2: Filter 5' end
    echo -e "\${SRA}_R1.fastq.bam\\n\${SRA}_R2.fastq.bam" | parallel -j 2 --eta \
    "samtools view -h \$RAW_DIR/{} | perl \$FILTER | samtools view -Sb - > \$FILT_DIR/{}"

    #echo "### Step 3A: Pair reads & mapping quality filter"
    perl \$COMBINER \$FILT_DIR/\${SRA}_R1.fastq.bam \$FILT_DIR/\${SRA}_R2.fastq.bam samtools \$MAPQ_FILTER | samtools view -bS -t \$FAIDX - | samtools sort -@ \$CPU -o \$TMP_DIR/\$SRA.bam -

    #echo "### Step 3.B: Add read group"
    java -jar /g/data/if89/apps/picard/2.27.4/picard.jar AddOrReplaceReadGroups INPUT=\$TMP_DIR/\$SRA.bam OUTPUT=\$PAIR_DIR/\$SRA.bam ID=\$SRA LB=\$SRA SM=\$LABEL PL=ILLUMINA PU=none

    #echo "### Step 4: Mark duplicates"
    java -jar /g/data/if89/apps/picard/2.27.4/picard.jar MarkDuplicates INPUT=\$PAIR_DIR/\$SRA.bam OUTPUT=\$REP_DIR/\$REP_LABEL.bam METRICS_FILE=\$REP_DIR/metrics.\$REP_LABEL.txt TMP_DIR=\$TMP_DIR ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=TRUE

    samtools index \$REP_DIR/\$REP_LABEL.bam

    cd dedup 
    mv \${label}_rep1.bam \${label}_ArimaHiC.bam 
    mv \${label}_rep1.bam.bai \${label}_ArimaHiC.bam.bai
    mv \${label}_ArimaHiC.bam \${output}
    mv \${label}_ArimaHiC.bam.bai \${output}

    """
}

process yahs {

    publishDir "${outputdir}", mode: 'copy'

    input: 
    tuple val(sample), val(fasta), val(R1), val(R2), val(bam), val(outputdir)

    output:
    tuple val(sample), path("${sample}_YAHS.fasta"), val(R1), val(R2), val(bam), val(outputdir)

    script:
    """
    module load yahs samtools
    samtools faidx ${fasta}
    yahs -e GATC,GANTC,CTNAG,TTAA -r 10000,20000,50000,100000,200000,500000,1000000 -o ${sample} ${fasta} ${bam}
    mv ${sample}_scaffolds_final.fa ${sample}_YAHS.fasta

    """
}

process generate_hicmap {

    publishDir "${outputdir}", mode: 'copy'

    input:
    tuple val(sample), path(scaffolds), val(R1), val(R2), val(bam), val(outputdir)

    output:
    tuple path("*HifiASM_YAHS*")



    script:
    """

    REF=${scaffolds}
    R1=${R1}
    R2=${R2}
    output=\${PWD}


    module load seqkit bwa python3 parallel
    source /g/data/xl04/ka6418/miniconda/etc/profile.d/conda.sh
    conda activate autohic

    label=${sample}
    WORK_DIR=\${PBS_JOBFS}
    mkdir -p \${WORK_DIR}/fastq
    ln -s \${R1} \${WORK_DIR}/fastq/\${label}_R1.fastq.gz
    ln -s \${R2} \${WORK_DIR}/fastq/\${label}_R2.fastq.gz

  
    cd \${WORK_DIR}
    seqkit sort -lr \${REF} > base.fa
    bwa index base.fa > bwa_index.log 
    wait
    seqkit fx2tab -nl base.fa > base.sizes
    JUICER_DIR=/g/data/xl04/ka6418/HiC_Assembly_Experiement/juicer
    python \$JUICER_DIR/misc/generate_site_positions.py Arima base \${WORK_DIR}/base.fa
    wait
    echo "Generate site positions done"
    path_juicer=/g/data/xl04/ka6418/HiC_Assembly_Experiement/juicer
    path_3d=/g/data/xl04/ka6418/HiC_Experiment/tools/3d-dna
    path_greenhill=/g/data/xl04/ka6418/greenhill/GreenHill
    echo "Starting Juicer"
    bash /g/data/xl04/ka6418/HiC_Assembly_Experiement/juicer/scripts/juicer.sh -D \$path_juicer -d \${WORK_DIR} -y \${WORK_DIR}/base_Arima.txt -t \${PBS_NCPUS}  -g base -s Arima -z base.fa -p base.sizes >juicer.log.o 2>juicer.log.e
    wait
    echo "Juicer finished"
    echo "Running 3D-DNA"
    awk -f \$path_3d/utils/generate-assembly-file-from-fasta.awk base.fa >base.assembly 2>generate.log.e
    echo "finisheed .assembly"
    wait
    bash \$path_3d/visualize/run-assembly-visualizer.sh base.assembly aligned/merged_nodups.txt || true
    echo "finisheed .hic"
    wait
    python \$path_greenhill/utils/fasta_to_juicebox_assembly.py base.fa >base.ctg_info.assembly
    echo "finisheed final"
    echo "moving things now"
    mv base.hic \${label}_HifiASM_YAHS.hic
    mv base.ctg_info.assembly \${label}_HifiASM_YAHS.ctg_info.assembly
    mv \${label}_HifiASM_YAHS.hic \${output}
    mv \${label}_HifiASM_YAHS.ctg_info.assembly \${output}
    echo "finished move"
    """

}

workflow {
    inputCh.view()
    arimaCh = arima(inputCh)
    yahsCh = yahs(arimaCh)
    generate_hicmapCh = generate_hicmap(yahsCh)


}
