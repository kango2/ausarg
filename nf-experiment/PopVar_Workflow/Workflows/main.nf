/**
 * This Nextflow script describes a whole-genome analysis workflow for human genomic data.
 * The workflow carries out the following primary steps:
 * 1. Aligns sequencing reads to a reference genome using BWA.
 * 2. Marks and sorts duplicate reads.
 * 3. Base quality recalibration.
 * 4. Calls variants using the HaplotypeCaller.
 * 5. Aggregates individual GVCFs into a combined GVCF.
 * 6. Calls genotypes on the combined GVCF.
 * 7. Variant recalibration for SNPs and INDELs.
 * 8. Variant annotation for high-quality variants.
 * 
 * Input data and resources are set using the 'params' object, which allows the user to
 * easily configure paths to reference files, resource directories, and other essential parameters.
 * Furthermore, specific processes are included from module files to keep the main workflow 
 * organized and modular.
 * 
 * The workflow provides a summary of its execution once completed, including the duration and status.
 * 
 * It's recommended to run this workflow on a high-performance computing environment due to 
 * the computational demands of processing whole-genome sequencing data.
 
 * Authors: Kosar Hooshmand
 */
 
import groovy.json.JsonBuilder

// Define the pipeline input parameters (with default values)
params.male_fasta = "/g/data/te53/humanreference/GRCh38.p13/refresource/GRCh38.p13.alt.decoyhs38d1.hla.phix.sequin.patches.xy.fa"
params.female_fasta = "/g/data/te53/humanreference/GRCh38.p13/refresource/GRCh38.p13.alt.decoyhs38d1.hla.phix.sequin.patches.xo.fa"
params.male_alt = "/g/data/te53/humanreference/GRCh38.p13/refresource/GRCh38.p13.alt.decoyhs38d1.hla.phix.sequin.patches.xy.fa"
params.female_alt = "/g/data/te53/humanreference/GRCh38.p13/refresource/GRCh38.p13.alt.decoyhs38d1.hla.phix.sequin.patches.xo.fa"
params.male_fai = "/g/data/te53/humanreference/GRCh38.p13/refresource/GRCh38.p13.alt.decoyhs38d1.hla.phix.sequin.patches.xy.fa.fai"
params.female_fai = "/g/data/te53/humanreference/GRCh38.p13/refresource/GRCh38.p13.alt.decoyhs38d1.hla.phix.sequin.patches.xo.fa.fai"
params.male_dict = "/g/data/te53/humanreference/GRCh38.p13/refresource/GRCh38.p13.alt.decoyhs38d1.hla.phix.sequin.patches.xy.dict"
params.female_dict = "/g/data/te53/humanreference/GRCh38.p13/refresource/GRCh38.p13.alt.decoyhs38d1.hla.phix.sequin.patches.xo.dict"
params.malebaseName = "genomexy"
params.femalebaseName = "genomexo"
params.postaltjs = "/g/data/te53/software/bwa/0.7.17-r1198-dirty/bwa-postalt.js"
params.pl = "illumina"
params.pm = "hiseq"
params.bwaindexDir = "/g/data/te53/kh3349/Nextflow_project/BWA_index"
params.baseDirPath = "/g/data/te53/kh3349/Nextflow_project/Test_project/NCIG_Population_Project/GIAB-hg004"
params.dbSNP = "/g/data/te53/referencedata/openaccess/dbsnp154/GCF_000001405.38.seqidmod.vcf.gz"
params.hapmap = "/g/data/te53/kh3349/Nextflow_project/Resource/hapmap_3.3.hg38.vcf.gz"
params.omni = "/g/data/te53/kh3349/Nextflow_project/Resource/1000G_omni2.5.hg38.vcf.gz"
params.thousandG = "/g/data/te53/kh3349/Nextflow_project/Resource/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
params.mills = "/g/data/te53/kh3349/Nextflow_project/Resource/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
params.intervalDir = "/g/data/te53/kh3349/Nextflow_project/Intervals"
params.hapmap_index = "/g/data/te53/kh3349/Nextflow_project/Resource/hapmap_3.3.hg38.vcf.gz.tbi"
params.omni_index = "/g/data/te53/kh3349/Nextflow_project/Resource/1000G_omni2.5.hg38.vcf.gz.tbi"
params.thousandG_index = "/g/data/te53/kh3349/Nextflow_project/Resource/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi"
params.tbi = "/g/data/te53/referencedata/openaccess/dbsnp154/GCF_000001405.38.seqidmod.vcf.gz.tbi"
params.mills_index = "/g/data/te53/kh3349/Nextflow_project/Resource/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi"
params.minmem = 16
params.mem = 32
params.threads = 8
params.outputDir ="/g/data/te53/kh3349/Nextflow_project/Test_project/NCIG_Population_Project/GIAB-hg004/nextflow-outputs"

// Include specific process module files

include { 
    bwa_mem;
    sort_markdup;
    baserecalibrator;
    applyBQSR;
    haplotypecaller;
    makelist;
    combineGVCFs_male;
    combineGVCFs_female;
    genotypeGVCFs_male;
    genotypeGVCFs_female;
    variantrecalibrator_SNP_male;
    variantrecalibrator_SNP_female;
    variantrecalibrator_INDEL_male;
    variantrecalibrator_INDEL_female;
    applyVQSR_SNP_male;
    applyVQSR_SNP_female;
    applyVQSR_INDEL_male;
    applyVQSR_INDEL_female;
    variantannotator_SNP_male;
    variantannotator_SNP_female;
    variantannotator_INDEL_male;
    variantannotator_INDEL_female
} from "/g/data/te53/kh3349/Nextflow_project/Test_project/NCIG_Population_Project/GIAB-hg004/modules/local/wg-human-snp.nf"

include { 
    clean_work_files
} from "/g/data/te53/kh3349/Nextflow_project/Test_project/NCIG_Population_Project/GIAB-hg004/modules/local/utilities.nf"

include { 
    fetchdataFromDB
} from "/g/data/te53/kh3349/Nextflow_project/Test_project/NCIG_Population_Project/GIAB-hg004/modules/local/fetchFromDB.nf"

include { 
    makesexlists
} from "/g/data/te53/kh3349/Nextflow_project/Test_project/NCIG_Population_Project/GIAB-hg004/modules/local/make-sex-lists.nf"


// The main workflow starts here
Channel.value("dummy_value").set { trigger }

workflow {

    // Run the fetchdataFromDB process using the trigger
    fetchdataFromDB(trigger)

    // Process the output
    fetchdataFromDB.out
    .splitCsv(sep: '\t', header: false)
    .map { fastq1, fastq2, sex ->
        // Extract the desired ID from the file name using a regex pattern
        def matcher = (fastq1 =~ /^.*\/giab\/HG001\/(.*)(?:_R[12]\.fastq\.gz)$/)
        def fastq1BaseName = matcher ? matcher[0][1] : ""

        // Depending on the gender, select different reference files
        def alt = (sex == 'Male') ? params.male_alt : params.female_alt
        def baseName = (sex == 'Male') ? params.malebaseName : params.femalebaseName
        return tuple(fastq1BaseName, file(fastq1), file(fastq2), file(alt), file(params.postaltjs), params.pl, params.pm, params.bwaindexDir, baseName)
    }
    .set { bwa_inputs }

    // Validate the existence of the files before proceeding to the next step
    bwa_inputs.view().map { id, fastq1, fastq2, genome_alt, postalt, pl, pm, bwa_index, genomeBaseName ->
        assert fastq1.exists() : "File $fastq1 does not exist"
        assert fastq2.exists() : "File $fastq2 does not exist"
        assert genome_alt.exists() : "File $genome_alt does not exist"
        assert postalt.exists() : "File $postalt does not exist"
        return tuple(id, fastq1, fastq2, genome_alt, postalt, pl, pm, bwa_index, genomeBaseName)
    }

    // Call the bwa_mem process
    bwa_mem(bwa_inputs)

    // Split the fetched data again and extract the necessary information for sorting and marking duplicates
    fetchdataFromDB.out
    .splitCsv(sep: '\t', header: false)
    .map { fastq1, _, sex ->
        def matcher = (fastq1 =~ /^.*\/giab\/HG001\/(.*)(?:_R[12]\.fastq\.gz)$/)
        def fastq1BaseName = matcher ? matcher[0][1] : ""
        def fasta = (sex == 'Male') ? params.male_fasta : params.female_fasta
        tuple(fastq1BaseName, file(fasta))
    }
    .set { sort_markdup_fasta_channel }

    // Join the results of bwa-mem alignment with the fasta files
    bwa_mem.out
    .join(sort_markdup_fasta_channel, by: [0])
    .map { id, sam, fasta -> tuple(id, sam, fasta) }
    .set { sort_markdup_inputs }

    // Add validation for the inputs of the sort_markdup process
    sort_markdup_inputs.view().map { id, sam, fasta ->
        assert sam.exists() : "File $sam does not exist"
        assert fasta.exists() : "File $fasta does not exist"
        return tuple(id, sam, fasta)
    }
    
    // Call the sort_markdup process
    sort_markdup(sort_markdup_inputs)

    // After marking duplicates, join the results and select files with .sam extension for cleanup
    bwa_mem.out
    .join(sort_markdup.out, by: [0])
    .flatten()
    .filter{ it =~ /.sam$/ }
    .set{ cleanable_sams }
    clean_work_files(
        cleanable_sams)

    // Split the output of fetchdataFromDB on tabs and process each row to extract necessary information
    fetchdataFromDB.out
    .splitCsv(sep: '\t', header: false)
    .map { fastq1, _, sex ->
    def matcher = (fastq1 =~ /^.*\/giab\/HG001\/(.*)(?:_R[12]\.fastq\.gz)$/)
    def fastq1BaseName = matcher ? matcher[0][1] : ""
   
    def fasta = (sex == 'Male') ? params.male_fasta : params.female_fasta
    def fai = (sex == 'Male') ? params.male_fai : params.female_fai
    def dict = (sex == 'Male') ? params.male_dict : params.female_dict
    tuple(fastq1BaseName, sex, file(fasta), file(fai), file(dict))
    }
    .set { processed_data }
    
    // Map the processed data to include dbSNP and its index for base recalibration
    processed_data.map { id, sex, fasta, fai, dict ->
    tuple(id, fasta, fai, dict, file(params.dbSNP), file(params.tbi))
    }
    .set { baserecalibrator_fasta_channel }
    
    // Join the results of sort_markdup alignment with the reference files for base recalibration
    sort_markdup.out
    .join(baserecalibrator_fasta_channel, by: [0])
    .map { id, bam, bai, fasta, fai, dict, dbSNP, tbi -> 
    tuple(id, bam, bai, fasta, fai, dict, dbSNP, tbi) 
    }
    .set { baserecalibrator_inputs }
    
    // Call the baserecalibrator process
    baserecalibrator(baserecalibrator_inputs)

    // Join the outputs of sort_markdup and baserecalibrator for subsequent steps
    joinedCh = sort_markdup.out.join(baserecalibrator.out)

    // Map the processed data for inputs to applyBQSR
    processed_data.map { id, sex, fasta, fai, dict ->
    tuple(id, fasta, fai, dict)
    }
    .set { applyBQSR_fasta_channel }
    
    // Join the channels to prepare inputs for the ApplyBQSR step
    joinedCh
    .join(applyBQSR_fasta_channel, by: [0])
    .map { id, bam, bai, data, fasta, fai, dict -> 
    tuple(id, bam, bai, data, fasta, fai, dict) 
    }
    .set { applyBQSR_inputs }
    
    // Call the applyBQSR process
    applyBQSR(applyBQSR_inputs)
    
    // Map the processed data for inputs to HaplotypeCaller
    processed_data.map { id, sex, fasta, fai, dict ->
    tuple(id, sex, fasta, fai, dict, file(params.intervalDir))
    }
    .set { haplotypecaller_fasta_channel }

   // Join the channels to prepare inputs for the HaplotypeCaller step
   applyBQSR.out
   .join(haplotypecaller_fasta_channel, by: [0])
   .map { id, input_bam, input_bam_bai, sex, fasta, fai, dict, interval ->
    tuple(id, input_bam, input_bam_bai, sex, fasta, fai, dict, interval)
    }
    .set { haplotypecaller_inputs }

    // Call the haplotypecaller process
    haplotypecaller(haplotypecaller_inputs)
    // Call the makelist process
    makelist(haplotypecaller.out[0].collect())
    // Call the makesexlists process
    makesexlists(makelist.out)

    // Filter the output channel to exclude empty male lists
    filtered_male_list = makesexlists.out.male.filter { file -> file.size() > 0 }
    filtered_female_list = makesexlists.out.female.filter { file -> file.size() > 0 }

    // Conditionally run combineGVCFs_female and combineGVCFs_male
    if (filtered_male_list) {
        combineGVCFs_male( filtered_male_list, params.male_fasta, params.male_fai, params.male_dict)
    }
    if (filtered_female_list) {
    combineGVCFs_female( filtered_female_list, params.female_fasta, params.female_fai, params.female_dict)
    }
    
    // Filter outputs of combineGVCFs_female and combineGVCFs_male
    filtered_combineGVCFs_male_out = combineGVCFs_male.out[0].filter { file -> file.size() > 0 }
    filtered_combineGVCFs_female_out = combineGVCFs_female.out[0].filter { file -> file.size() > 0 }
    
    // Conditionally run genotypeGVCFs_female and genotypeGVCFs_male
    if (filtered_combineGVCFs_male_out) {
    genotypeGVCFs_male(filtered_combineGVCFs_male_out, combineGVCFs_male.out[1], params.male_fasta, params.male_fai, params.male_dict)
    }
    if (filtered_combineGVCFs_female_out) {
    genotypeGVCFs_female(filtered_combineGVCFs_female_out, combineGVCFs_female.out[1], params.female_fasta, params.female_fai, params.female_dict)
    }

    // Filter outputs of genotypeGVCFs_female and genotypeGVCFs_male
    filtered_genotypeGVCFs_male_out = genotypeGVCFs_male.out[0].filter { file -> file.size() > 0 }
    filtered_genotypeGVCFs_female_out = genotypeGVCFs_female.out[0].filter { file -> file.size() > 0 }
    
    // Conditionally run variantrecalibrator_SNP_female and variantrecalibrator_SNP_male
    if (filtered_genotypeGVCFs_male_out) {
    variantrecalibrator_SNP_male(filtered_genotypeGVCFs_male_out, genotypeGVCFs_male.out[1], params.male_fasta, params.male_fai, params.male_dict, params.dbSNP, params.tbi, params.hapmap, params.hapmap_index, params.omni, params.omni_index, params.thousandG, params.thousandG_index)
    }
    if (filtered_genotypeGVCFs_female_out) {
    variantrecalibrator_SNP_female(filtered_genotypeGVCFs_female_out, genotypeGVCFs_female.out[1], params.female_fasta, params.female_fai, params.female_dict, params.dbSNP, params.tbi, params.hapmap, params.hapmap_index, params.omni, params.omni_index, params.thousandG, params.thousandG_index)
    }
    
    // Conditionally run variantrecalibrator_INDEL_female and variantrecalibrator_INDEL_male
    if (filtered_genotypeGVCFs_male_out) {
    variantrecalibrator_INDEL_male(filtered_genotypeGVCFs_male_out, genotypeGVCFs_male.out[1], params.male_fasta, params.male_fai, params.male_dict, params.dbSNP, params.tbi, params.mills, params.mills_index)
    }
    if (filtered_genotypeGVCFs_female_out) {
    variantrecalibrator_INDEL_female(filtered_genotypeGVCFs_female_out, genotypeGVCFs_female.out[1], params.female_fasta, params.female_fai, params.female_dict, params.dbSNP, params.tbi, params.mills, params.mills_index)
    }

    // Filter outputs of variantrecalibrator_SNP_female and variantrecalibrator_SNP_male
    filtered_variantrecalibrator_SNP_male_out = variantrecalibrator_SNP_male.out[0].filter { file -> file.size() > 0 }
    filtered_variantrecalibrator_SNP_female_out = variantrecalibrator_SNP_female.out[0].filter { file -> file.size() > 0 }

    // Conditionally run applyVQSR_SNP_female and applyVQSR_SNP_male
    if (filtered_variantrecalibrator_SNP_male_out) {
    applyVQSR_SNP_male(filtered_genotypeGVCFs_male_out, genotypeGVCFs_male.out[1], filtered_variantrecalibrator_SNP_male_out, variantrecalibrator_SNP_male.out[1], variantrecalibrator_SNP_male.out[2], params.male_fasta, params.male_fai, params.male_dict)
    }
    if (filtered_variantrecalibrator_SNP_female_out) {
    applyVQSR_SNP_female(filtered_genotypeGVCFs_female_out, genotypeGVCFs_female.out[1], filtered_variantrecalibrator_SNP_female_out, variantrecalibrator_SNP_female.out[1], variantrecalibrator_SNP_female.out[2], params.female_fasta, params.female_fai, params.female_dict)
    }

    // Filter outputs of variantrecalibrator_INDEL_female and variantrecalibrator_INDEL_male
    filtered_variantrecalibrator_INDEL_male_out = variantrecalibrator_INDEL_male.out[0].filter { file -> file.size() > 0 }
    filtered_variantrecalibrator_INDEL_female_out = variantrecalibrator_INDEL_female.out[0].filter { file -> file.size() > 0 }

    // Conditionally run applyVQSR_INDEL_female and applyVQSR_INDEL_male
    if (filtered_variantrecalibrator_INDEL_male_out) {
    applyVQSR_INDEL_male(filtered_genotypeGVCFs_male_out, genotypeGVCFs_male.out[1], filtered_variantrecalibrator_INDEL_male_out, variantrecalibrator_INDEL_male.out[1], variantrecalibrator_INDEL_male.out[2], params.male_fasta, params.male_fai, params.male_dict)
    }
    if (filtered_variantrecalibrator_INDEL_female_out) {
    applyVQSR_INDEL_female(filtered_genotypeGVCFs_female_out, genotypeGVCFs_female.out[1], filtered_variantrecalibrator_INDEL_female_out, variantrecalibrator_INDEL_female.out[1], variantrecalibrator_INDEL_female.out[2], params.female_fasta, params.female_fai, params.female_dict)
    }

    // Filter outputs of applyVQSR_SNP_female and applyVQSR_SNP_male
    filtered_applyVQSR_SNP_male_out = applyVQSR_SNP_male.out[0].filter { file -> file.size() > 0 }
    filtered_applyVQSR_SNP_female_out = applyVQSR_SNP_female.out[0].filter { file -> file.size() > 0 }

    // Conditionally run variantannotator_SNP_male and variantannotator_SNP_female
    if (filtered_applyVQSR_SNP_male_out) {
    variantannotator_SNP_male(filtered_applyVQSR_SNP_male_out, applyVQSR_SNP_male.out[1], params.dbSNP, params.tbi)
    }
    if (filtered_applyVQSR_SNP_female_out) {
    variantannotator_SNP_female(filtered_applyVQSR_SNP_female_out, applyVQSR_SNP_female.out[1], params.dbSNP, params.tbi)
    }

    // Filter outputs of applyVQSR_INDEL_female and applyVQSR_INDEL_male
    filtered_applyVQSR_INDEL_male_out = applyVQSR_INDEL_male.out[0].filter { file -> file.size() > 0 }
    filtered_applyVQSR_INDEL_female_out = applyVQSR_INDEL_female.out[0].filter { file -> file.size() > 0 }

    // Conditionally run variantannotator_INDEL_male and variantannotator_INDEL_female
    if (filtered_applyVQSR_INDEL_male_out) {
    variantannotator_INDEL_male(filtered_applyVQSR_INDEL_male_out, applyVQSR_INDEL_male.out[1], params.dbSNP, params.tbi)
    }
    if (filtered_applyVQSR_INDEL_female_out) {
    variantannotator_INDEL_female(filtered_applyVQSR_INDEL_female_out, applyVQSR_INDEL_female.out[1], params.dbSNP, params.tbi)
    }
}
