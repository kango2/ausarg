/*
* This script contains a series of processes to perform variant calling on genomic data. Here's an overview:

* 1. **bwa_mem**: Aligns reads using the BWA-MEM algorithm.
* 2. **sort_markdup**: Fixes mate information and sorts BAM/CRAM files, then marks duplicates.
* 3. **baserecalibrator**: Executes the GATK BaseRecalibrator tool to generate recalibration tables based on various known datasets.
* 4. **applyBQSR**: Applies the base quality score recalibration (BQSR) on the sorted BAM file.
* 5. **haplotypecaller**: Calls germline SNPs and indels via local re-assembly of haplotypes using the GATK HaplotypeCaller tool.
* 6. **makelist**: Generates a list of GVCF files.
* 7. **combineGVCFs**: Combines multiple GVCF files into one using the GATK CombineGVCFs tool.
* 8. **genotypeGVCFs**: Performs joint genotyping on GVCF files using the GATK GenotypeGVCFs tool.
* 9. **variantrecalibrator_SNP**: Executes the GATK VariantRecalibrator tool for SNP variant quality score recalibration (VQSR).
* 10. **applyVQSR_SNP**: Applies VQSR on SNP calls using the GATK ApplyVQSR tool.
* 11. **variantannotator_SNP**: Annotates SNP variants using the bcftools annotate tool.
* 12. **variantrecalibrator_INDEL**: Similar to the SNP variant recalibration, but this one focuses on INDELs.
* 13. **applyVQSR_INDEL**: Similar to the SNP Applies VQSR, but this one focuses on INDELs.
* 14. **variantannotator_INDEL**: Similar to the SNP VariantAnnotator, but this one focuses on INDELs.

* Each process uses specific modules and has specified inputs and outputs. Note that this script seems to focus on human genomics, using tools such as BWA for * alignment and GATK for variant calling and post-processing. Some of these processes also use reference genome files and known datasets for improved accuracy.

* Authors: Kosar Hooshmand
*/

import groovy.json.JsonBuilder

process bwa_mem {
    module 'bwa/0.7.17'
    module 'k8/0.2.5'
    
    label 'bwa_mem'
    errorStrategy 'ignore'

    input:
    tuple val(pair_id), path(fastq1), path(fastq2), file(genome_alt), file(postalt), val(pl), val(pm), path(bwa_index), val(genomeBaseName)

    output:
    tuple val(pair_id), path("${pair_id}.sam")

    script:
    """
    bwa mem -t 16 -Y -K 100000000 -R "@RG\\tID:${pair_id}\\tPL:${pl}\\tPM:${pm}\\tSM:${pair_id}" "${bwa_index}/${genomeBaseName}" ${fastq1} ${fastq2} | k8 ${postalt} -a ${genome_alt} > ${pair_id}.sam
    """
}

process sort_markdup {
    module 'samtools/1.12'

    publishDir "${params.outputDir}/${task.process}", mode: 'copy', overwrite: false
    label 'sort_markdup'
    errorStrategy 'ignore'
   
    input:
    tuple val(pair_id), path(sam), file(genome)

    output:
    tuple val(pair_id), path("${pair_id}.sortdedup.bam"), path("${pair_id}.sortdedup.bam.bai")

    script:
    def tmpdir = '/g/data/te53/kh3349/Nextflow_project/tmpdir'

    """
    samtools fixmate -m --threads 16 --reference ${genome} ${sam} - | \
    samtools sort --threads 16 --reference ${genome} -T ${tmpdir}/${pair_id}.sorttmp | \
    samtools markdup -T ${tmpdir}/${pair_id}.deduptmp -O BAM --threads 16 - ${pair_id}.sortdedup.bam
    samtools index ${pair_id}.sortdedup.bam
    """
}

process baserecalibrator {
    module 'java'
    module 'gatk/4.2.5.0'

    publishDir "${params.outputDir}/${task.process}", mode: 'copy', overwrite: false
    label 'baserecalibrator'
    errorStrategy 'ignore'

    input:
    tuple val(pair_id), path(bam), path(bam_bai), file(genome), file(genome_fai), file(genome_dict), file(dbSNP), file(dbSNP_index)

    output:
    tuple val(pair_id), path("${pair_id}.recal_data.table") 

    script: 

    """
    gatk --java-options '-XX:+UseSerialGC -Xms${params.minmem}g -Xmx${params.mem}g' BaseRecalibrator \
    -I ${bam} \
    -R ${genome} \
    --known-sites ${dbSNP} \
    -O ${pair_id}.recal_data.table 
    """
}

process applyBQSR {
    module 'java'
    module 'gatk/4.2.5.0'

    publishDir "${params.outputDir}/${task.process}", mode: 'copy', overwrite: false
    label 'applyBQSR'
    errorStrategy 'ignore'

    input:
    tuple val(pair_id), path(bam), path(bam_bai), path(data), file(genome), file(genome_fai), file(genome_dict)

    output:
    tuple val(pair_id), path("${pair_id}.recalibrated.bam"), path("${pair_id}.recalibrated.bai")

    script: 

    """
    gatk --java-options '-XX:+UseSerialGC -Xms${params.minmem}g -Xmx${params.mem}g' ApplyBQSR \
    -R ${genome} \
    -I ${bam} \
    --bqsr-recal-file ${data} \
    -O ${pair_id}.recalibrated.bam
    """
}

process haplotypecaller {
    module 'java'
    module 'gatk/4.2.5.0'

    publishDir "${params.outputDir}/${task.process}", mode: 'copy', overwrite: false
    label 'haplotypecaller'
    errorStrategy 'ignore'

    input:
    tuple val(pair_id), path(input_bam), path(input_bam_bai), val(sex), file(genome), file(genome_fai), file(genome_dict), path(interval_files)
    
    output:
    tuple val(pair_id), path("${pair_id}_${sex}.raw_variants.g.vcf.gz")
    tuple val(pair_id), path("${pair_id}_${sex}.raw_variants.g.vcf.gz.tbi")

    script:

    """
    gatk --java-options '-XX:+UseSerialGC -Xms${params.minmem}g -Xmx${params.mem}g' HaplotypeCaller \
    -R ${genome} \
    -I ${input_bam} \
    --native-pair-hmm-threads ${params.threads} \
    --sample-ploidy 2 \
    --min-base-quality-score 20 \
    --gvcf-gq-bands 10 \
    --gvcf-gq-bands 20 \
    --gvcf-gq-bands 30 \
    --gvcf-gq-bands 60 \
    -O "${pair_id}_${sex}.raw_variants.g.vcf.gz" \
    -ERC GVCF \
    -L "${interval_files}/0000-scattered.interval_list" -L "${interval_files}/0001-scattered.interval_list" -L "${interval_files}/0002-scattered.interval_list" -L "${interval_files}/0003-scattered.interval_list" -L "${interval_files}/0004-scattered.interval_list" -L "${interval_files}/0005-scattered.interval_list" -L "${interval_files}/0006-scattered.interval_list" -L "${interval_files}/0007-scattered.interval_list" -L "${interval_files}/0008-scattered.interval_list" -L "${interval_files}/0009-scattered.interval_list"
    """
}

process makelist {

    publishDir "${params.outputDir}/${task.process}", mode: 'copy', overwrite: false
    label 'makelist'
    errorStrategy 'ignore'

    input:
    val gvcf_file

    output:
    path("gvcf.list")

    script:
    
    """
    echo "${gvcf_file.join('\n')}" | awk '!(NR%2)' > gvcf.list
    """
}

process combineGVCFs_male {
    module 'java'
    module 'gatk/4.2.5.0'

    publishDir "${params.outputDir}/${task.process}", mode: 'copy', overwrite: false
    label 'combineGVCFs_male'
    errorStrategy 'ignore'

    input:
    file male_list 
    file params.male_fasta
    file params.male_fai
    file params.male_dict

    output:
    path("combined_male.g.vcf.gz") 
    path("combined_male.g.vcf.gz.tbi")

    script:
    def tmpdir = '/g/data/te53/kh3349/Nextflow_project/tmpdir'

    """
    gatk --java-options '-XX:+UseSerialGC -Xms${params.minmem}g -Xmx${params.mem}g' CombineGVCFs \
    -R ${params.male_fasta} \
    -V ${male_list} \
    -O combined_male.g.vcf.gz \
    -G StandardAnnotation \
    -G AS_StandardAnnotation \
    --tmp-dir ${tmpdir} \
    """
}

process combineGVCFs_female {
    module 'java'
    module 'gatk/4.2.5.0'

    publishDir "${params.outputDir}/${task.process}", mode: 'copy', overwrite: false
    label 'combineGVCFs_female'
    errorStrategy 'ignore'

    input:
    file female_list 
    file params.female_fasta
    file params.female_fai
    file params.female_dict

    output:
    path("combined_female.g.vcf.gz") 
    path("combined_female.g.vcf.gz.tbi")

    script:
    def tmpdir = '/g/data/te53/kh3349/Nextflow_project/tmpdir'

    """
    gatk --java-options '-XX:+UseSerialGC -Xms${params.minmem}g -Xmx${params.mem}g' CombineGVCFs \
    -R ${params.female_fasta} \
    -V ${female_list} \
    -O combined_female.g.vcf.gz \
    -G StandardAnnotation \
    -G AS_StandardAnnotation \
    --tmp-dir ${tmpdir} \
    """
}

process genotypeGVCFs_male {
    module 'java'
    module 'gatk/4.2.5.0'

    publishDir "${params.outputDir}/${task.process}", mode: 'copy', overwrite: false
    label 'genotypeGVCFs_male'
    errorStrategy 'ignore'

    input:
    file combined_male_GVCF
    file combined_male_GVCF_index
    file params.male_fasta
    file params.male_fai
    file params.male_dict

    output:
    path("combined_male_raw.vcf.gz")
    path("combined_male_raw.vcf.gz.tbi")
    
    script:
    def tmpdir = '/g/data/te53/kh3349/Nextflow_project/tmpdir'

    """
    gatk --java-options '-XX:+UseSerialGC -Xms${params.minmem}g -Xmx${params.mem}g' GenotypeGVCFs \
    -R ${params.male_fasta} \
    -V ${combined_male_GVCF} \
    -stand-call-conf 50 \
    -A Coverage \
    -A FisherStrand \
    -A StrandOddsRatio \
    -A MappingQualityRankSumTest \
    -A QualByDepth \
    -A RMSMappingQuality \
    -A ReadPosRankSumTest \
    --allow-old-rms-mapping-quality-annotation-data \
    -O combined_male_raw.vcf.gz \
    --tmp-dir ${tmpdir}
    """
}

process genotypeGVCFs_female {
    module 'java'
    module 'gatk/4.2.5.0'

    publishDir "${params.outputDir}/${task.process}", mode: 'copy', overwrite: false
    label 'genotypeGVCFs_female'
    errorStrategy 'ignore'

    input:
    file combined_female_GVCF
    file combined_female_GVCF_index
    file params.female_fasta
    file params.female_fai
    file params.female_dict

    output:
    path("combined_female_raw.vcf.gz")
    path("combined_female_raw.vcf.gz.tbi")
    
    script:
    def tmpdir = '/g/data/te53/kh3349/Nextflow_project/tmpdir'

    """
    gatk --java-options '-XX:+UseSerialGC -Xms${params.minmem}g -Xmx${params.mem}g' GenotypeGVCFs \
    -R ${params.female_fasta} \
    -V ${combined_female_GVCF} \
    -stand-call-conf 50 \
    -A Coverage \
    -A FisherStrand \
    -A StrandOddsRatio \
    -A MappingQualityRankSumTest \
    -A QualByDepth \
    -A RMSMappingQuality \
    -A ReadPosRankSumTest \
    --allow-old-rms-mapping-quality-annotation-data \
    -O combined_female_raw.vcf.gz \
    --tmp-dir ${tmpdir}
    """
}

process variantrecalibrator_SNP_male {
    module 'java'
    module 'gatk/4.2.5.0'
    
    publishDir "${params.outputDir}/${task.process}", mode: 'copy', overwrite: false
    label 'variantrecalibrator_SNP_male'
    errorStrategy 'ignore'
    
    input:
    file input_male_VCF
    file input_male_VCF_index
    file params.male_fasta
    file params.male_fai
    file params.male_dict
    file params.dbSNP
    file params.tbi
    file params.hapmap
    file params.hapmap_index
    file params.omni 
    file params.omni_index
    file params.thousandG
    file params.thousandG_index
    
    output:
    path("SNP_male.recal")
    path("SNP_male.recal.idx")
    path("SNP_male.tranches")
    
    script:
    
    """
    gatk --java-options '-XX:+UseSerialGC -Xms${params.minmem}g -Xmx${params.mem}g' VariantRecalibrator \
    -R ${params.male_fasta} \
    -V ${input_male_VCF} \
    -AS \
    -tranche 100.0 \
    -tranche 99.95 \
    -tranche 99.9 \
    -tranche 99.8 \
    -tranche 99.6 \
    -tranche 99.5 \
    -tranche 99.4 \
    -tranche 99.3 \
    -tranche 99.0 \
    -tranche 98.0 \
    -tranche 97.0 \
    -tranche 90.0 \
    -an QD -an ReadPosRankSum -an FS -an SOR -an DP \
    --mode SNP \
    --max-gaussians 1 \
    -resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${params.hapmap} \
    -resource:omni,known=false,training=true,truth=false,prior=12.0 ${params.omni} \
    -resource:1000G,known=false,training=true,truth=false,prior=10.0 ${params.thousandG} \
    -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${params.dbSNP} \
    -O SNP_male.recal \
    --tranches-file SNP_male.tranches 
  """
}

process variantrecalibrator_SNP_female {
    module 'java'
    module 'gatk/4.2.5.0'
    
    publishDir "${params.outputDir}/${task.process}", mode: 'copy', overwrite: false
    label 'variantrecalibrator_SNP_female'
    errorStrategy 'ignore'
    
    input:
    file input_female_VCF
    file input_female_VCF_index
    file params.female_fasta
    file params.female_fai
    file params.female_dict
    file params.dbSNP
    file params.tbi
    file params.hapmap
    file params.hapmap_index
    file params.omni 
    file params.omni_index
    file params.thousandG
    file params.thousandG_index
    
    output:
    path("SNP_female.recal")
    path("SNP_female.recal.idx")
    path("SNP_female.tranches")
    
    script:
    
    """
    gatk --java-options '-XX:+UseSerialGC -Xms${params.minmem}g -Xmx${params.mem}g' VariantRecalibrator \
    -R ${params.female_fasta} \
    -V ${input_female_VCF} \
    -AS \
    -tranche 100.0 \
    -tranche 99.95 \
    -tranche 99.9 \
    -tranche 99.8 \
    -tranche 99.6 \
    -tranche 99.5 \
    -tranche 99.4 \
    -tranche 99.3 \
    -tranche 99.0 \
    -tranche 98.0 \
    -tranche 97.0 \
    -tranche 90.0 \
    -an QD -an ReadPosRankSum -an FS -an SOR -an DP \
    --mode SNP \
    --max-gaussians 4 \
    -resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${params.hapmap} \
    -resource:omni,known=false,training=true,truth=false,prior=12.0 ${params.omni} \
    -resource:1000G,known=false,training=true,truth=false,prior=10.0 ${params.thousandG} \
    -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${params.dbSNP} \
    -O SNP_female.recal \
    --tranches-file SNP_female.tranches 
  """
}

process variantrecalibrator_INDEL_male {
    module 'java'
    module 'gatk/4.2.5.0'
    
    publishDir "${params.outputDir}/${task.process}", mode: 'copy', overwrite: false
    label 'variantrecalibrator_INDEL_male'
    errorStrategy 'ignore'
    
    input:
    file input_male_VCF
    file input_male_VCF_index
    file params.male_fasta
    file params.male_fai
    file params.male_dict
    file params.dbSNP
    file params.tbi
    file params.mills
    file params.mills_index
    
    output:
    path("INDEL_male.recal")
    path("INDEL_male.recal.idx") 
    path("INDEL_male.tranches") 
    
    script:
    
    """
    gatk --java-options '-XX:+UseSerialGC -Xms${params.minmem}g -Xmx${params.mem}g' VariantRecalibrator \
    -R ${params.male_fasta} \
    -V ${input_male_VCF} \
    -resource:mills,known=false,training=true,truth=true,prior=12.0 ${params.mills} \
    -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${params.dbSNP} \
    -AS \
    -an QD -an DP -an FS -an SOR -an ReadPosRankSum \
    --mode INDEL \
    -tranche 100.0 \
    -tranche 99.9 \
    -tranche 99.0 \
    -tranche 90.0 \
    --max-gaussians 1 \
    -O INDEL_male.recal \
    --tranches-file INDEL_male.tranches 
    """
}

process variantrecalibrator_INDEL_female {
    module 'java'
    module 'gatk/4.2.5.0'
    
    publishDir "${params.outputDir}/${task.process}", mode: 'copy', overwrite: false
    label 'variantrecalibrator_INDEL_female'
    errorStrategy 'ignore'
    
    input:
    file input_female_VCF
    file input_female_VCF_index
    file params.female_fasta
    file params.female_fai
    file params.female_dict
    file params.dbSNP
    file params.tbi
    file params.mills
    file params.mills_index
    
    output:
    path("INDEL_female.recal")
    path("INDEL_female.recal.idx") 
    path("INDEL_female.tranches") 
    
    script:
    
    """
    gatk --java-options '-XX:+UseSerialGC -Xms${params.minmem}g -Xmx${params.mem}g' VariantRecalibrator \
    -R ${params.female_fasta} \
    -V ${input_female_VCF} \
    -resource:mills,known=false,training=true,truth=true,prior=12.0 ${params.mills} \
    -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${params.dbSNP} \
    -AS \
    -an QD -an DP -an FS -an SOR -an ReadPosRankSum \
    --mode INDEL \
    -tranche 100.0 \
    -tranche 99.9 \
    -tranche 99.0 \
    -tranche 90.0 \
    --max-gaussians 1 \
    -O INDEL_female.recal \
    --tranches-file INDEL_female.tranches 
    """
}
process applyVQSR_SNP_male {
    module 'java'
    module 'gatk/4.2.5.0'

    publishDir "${params.outputDir}/${task.process}", mode: 'copy', overwrite: false
    label 'applyVQSR_SNP_male'
    errorStrategy 'ignore'

    input:
    file input_male_VCF
    file input_male_VCF_index
    file snp_male_recal
    file snp_male_recal_index
    file snp_male_tranches
    file params.male_fasta
    file params.male_fai
    file params.male_dict

    output:
    path("SNP_male_recal.vcf.gz") 
    path("SNP_male_recal.vcf.gz.tbi")
    
    script:

    """
    gatk --java-options '-XX:+UseSerialGC -Xms${params.minmem}g -Xmx${params.mem}g' ApplyVQSR \
    -R ${params.male_fasta} \
    --recal-file ${snp_male_recal} \
    --tranches-file ${snp_male_tranches} \
    --ts-filter-level 99.9 \
    -AS \
    -mode SNP \
    -V ${input_male_VCF} \
    -O SNP_male_recal.vcf.gz 
    """
}

process applyVQSR_SNP_female {
    module 'java'
    module 'gatk/4.2.5.0'

    publishDir "${params.outputDir}/${task.process}", mode: 'copy', overwrite: false
    label 'applyVQSR_SNP_female'
    errorStrategy 'ignore'

    input:
    file input_female_VCF
    file input_female_VCF_index
    file snp_female_recal
    file snp_female_recal_index
    file snp_female_tranches
    file params.female_fasta
    file params.female_fai
    file params.female_dict

    output:
    path("SNP_female_recal.vcf.gz") 
    path("SNP_female_recal.vcf.gz.tbi")
    
    script:

    """
    gatk --java-options '-XX:+UseSerialGC -Xms${params.minmem}g -Xmx${params.mem}g' ApplyVQSR \
    -R ${params.female_fasta} \
    --recal-file ${snp_female_recal} \
    --tranches-file ${snp_female_tranches} \
    --ts-filter-level 99.9 \
    -AS \
    -mode SNP \
    -V ${input_female_VCF} \
    -O SNP_female_recal.vcf.gz 
    """
}

process applyVQSR_INDEL_male {
    module 'java'
    module 'gatk/4.2.5.0'

    publishDir "${params.outputDir}/${task.process}", mode: 'copy', overwrite: false
    label 'applyVQSR_INDEL_male'
    errorStrategy 'ignore'

    input:
    file input_male_VCF
    file input_male_VCF_index
    file indel_male_recal
    file indel_male_recal_index
    file indel_male_tranches
    file params.male_fasta
    file params.male_fai
    file params.male_dict
    
    output:
    path("INDEL_male_recal.vcf.gz") 
    path("INDEL_male_recal.vcf.gz.tbi")
    
    script:
    
    """
    gatk --java-options '-XX:+UseSerialGC -Xms${params.minmem}g -Xmx${params.mem}g' ApplyVQSR \
    -R ${params.male_fasta} \
    --recal-file ${indel_male_recal} \
    --tranches-file ${indel_male_tranches} \
    --ts-filter-level 99.9 \
    -AS \
    -mode INDEL \
    -V ${input_male_VCF} \
    -O INDEL_male_recal.vcf.gz 
    """
}

process applyVQSR_INDEL_female {
    module 'java'
    module 'gatk/4.2.5.0'

    publishDir "${params.outputDir}/${task.process}", mode: 'copy', overwrite: false
    label 'applyVQSR_INDEL_female'
    errorStrategy 'ignore'

    input:
    file input_female_VCF
    file input_female_VCF_index
    file indel_female_recal
    file indel_female_recal_index
    file indel_female_tranches
    file params.female_fasta
    file params.female_fai
    file params.female_dict
    
    output:
    path("INDEL_female_recal.vcf.gz") 
    path("INDEL_female_recal.vcf.gz.tbi")
    
    script:
    
    """
    gatk --java-options '-XX:+UseSerialGC -Xms${params.minmem}g -Xmx${params.mem}g' ApplyVQSR \
    -R ${params.female_fasta} \
    --recal-file ${indel_female_recal} \
    --tranches-file ${indel_female_tranches} \
    --ts-filter-level 99.9 \
    -AS \
    -mode INDEL \
    -V ${input_female_VCF} \
    -O INDEL_female_recal.vcf.gz 
    """
}

process variantannotator_SNP_male {
    module 'bcftools'

    publishDir "${params.outputDir}/${task.process}", mode: 'copy', overwrite: false
    label 'variantannotator_SNP_male' 
    errorStrategy 'ignore'
    
    input:
    file vcf_male_SNP
    file vcf_male_SNP_index
    file params.dbSNP
    file params.tbi
    
    output:
    path("SNP_male_recal.annotated.vcf.gz")
    
    script:
    
    """
    bcftools annotate --columns ID  --annotations ${params.dbSNP} --threads ${params.threads} -Oz -o SNP_male_recal.annotated.vcf.gz ${vcf_male_SNP}
    """
}

process variantannotator_SNP_female {
    module 'bcftools'

    publishDir "${params.outputDir}/${task.process}", mode: 'copy', overwrite: false
    label 'variantannotator_SNP_female' 
    errorStrategy 'ignore'
    
    input:
    file vcf_female_SNP
    file vcf_female_SNP_index
    file params.dbSNP
    file params.tbi
    
    output:
    path("SNP_female_recal.annotated.vcf.gz")
    
    script:
    
    """
    bcftools annotate --columns ID  --annotations ${params.dbSNP} --threads ${params.threads} -Oz -o SNP_female_recal.annotated.vcf.gz ${vcf_female_SNP}
    """
}

process variantannotator_INDEL_male {
    module 'bcftools'

    publishDir "${params.outputDir}/${task.process}", mode: 'copy', overwrite: false
    label 'variantannotator_INDEL_male'
    errorStrategy 'ignore'
    
    input:
    file vcf_male_INDEL
    file vcf_male_INDEL_index
    file params.dbSNP
    file params.tbi
    
    output:
    path("INDEL_male_recal.annotated.vcf.gz")
    
    script:
    
    """
    bcftools annotate --columns ID  --annotations ${params.dbSNP} --threads ${params.threads} -Oz -o INDEL_male_recal.annotated.vcf.gz ${vcf_male_INDEL}
    """
}

process variantannotator_INDEL_female {
    module 'bcftools'

    publishDir "${params.outputDir}/${task.process}", mode: 'copy', overwrite: false
    label 'variantannotator_INDEL_female'
    errorStrategy 'ignore'
    
    input:
    file vcf_female_INDEL
    file vcf_female_INDEL_index
    file params.dbSNP
    file params.tbi
    
    output:
    path("INDEL_female_recal.annotated.vcf.gz")
    
    script:
    
    """
    bcftools annotate --columns ID  --annotations ${params.dbSNP} --threads ${params.threads} -Oz -o INDEL_female_recal.annotated.vcf.gz ${vcf_female_INDEL}
    """
}
