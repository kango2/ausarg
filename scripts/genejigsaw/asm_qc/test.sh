./asm_qc.sh \
-f /g/data/xl04/ka6418/bassiana/all_assemblies/rBasDup_HifiASM_YAHS_Optimised.fasta \
-i "/g/data/xl04/ka6418/LamDe/raw/illumina/dnaseq/354040_AusARG_UNSW_HJTCJDSX3_ATCGGCGAAG-TCCAAGAATT_S1_L004_R1_001.fastq.gz;/g/data/xl04/ka6418/LamDe/raw/illumina/dnaseq/354040_AusARG_UNSW_HJTCJDSX3_ATCGGCGAAG-TCCAAGAATT_S1_L004_R2_001.fastq.gz:/g/data/xl04/ka6418/LamDe/raw/illumina/dnaseq/354041_AusARG_UNSW_HJTCJDSX3_CCGTGACCGA-CCGAACGTTG_S2_L004_R1_001.fastq.gz;/g/data/xl04/ka6418/LamDe/raw/illumina/dnaseq/354041_AusARG_UNSW_HJTCJDSX3_CCGTGACCGA-CCGAACGTTG_S2_L004_R2_001.fastq.gz:/g/data/xl04/ka6418/LamDe/raw/illumina/dnaseq/354042_AusARG_UNSW_HJTCJDSX3_TATCCGAGGC-GGCCAATAAG_S3_L004_R1_001.fastq.gz;/g/data/xl04/ka6418/LamDe/raw/illumina/dnaseq/354042_AusARG_UNSW_HJTCJDSX3_TATCCGAGGC-GGCCAATAAG_S3_L004_R2_001.fastq.gz" \
-p /g/data/xl04/bpadata/Lampropholis_delicata/raw/pacbio/350750_AusARG_AGRF_DA064169.ccs.fq.gz:/g/data/xl04/bpadata/Lampropholis_delicata/raw/pacbio/350750_AusARG_AGRF_DA095606.ccs.fq.gz \
-o /g/data/xl04/ka6418/LamDe/raw/ont/lamde_ont.fastq.gz \
-d /g/data/xl04/ka6418/genejigsaw_testing/lamde_yahs \
-l /g/data/xl04/ka6418/genejigsaw_testing/lamde_yahs/logs \
-P xl04 \
-s gdata/if89+gdata/xl04
