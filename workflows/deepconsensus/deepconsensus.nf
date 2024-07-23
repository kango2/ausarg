//Output & Sample paths need to be a separate channel as well and then combine 


params.image = "/g/data/if89/shpcroot/containers/google/deepconsensus/1.2.0-gpu/google-deepconsensus-1.2.0-gpu-sha256:5a991780ecb7682bf7b93754a6774271aebedbd19b018b6a7aabc7f2162f8e88.sif"
params.chunks = 10 

process pbindex {

    input:

    val (subread)

    output:
    
    tuple val (subread) , path ("indexdone.txt")

    script:

    """
    module load singularity 
    singularity exec ${params.image} pbindex ${subread} -j \${PBS_NCPUS}
    touch "indexdone.txt"

    """
}


process ccs {

    input: 

    tuple val (subread) , path (pbi), val (id)

    output :
    
    tuple val (subread) , path ("${id}.ccs.bam") , val (id)

    script:

    """
    module load singularity 
    singularity exec ${params.image} ccs --min-rq=0.88 -j \${PBS_NCPUS} --chunk="${id}"/"${params.chunks}" ${subread} "${id}.ccs.bam"

    """

}


process actc 
{

    input:
    tuple val (subread) , path (ccs) , val  (id)

    output:
    tuple val (subread) , path ("${id}.subreads_to_ccs.bam") , path ("${ccs}") , val (id)

    script:
    """
    module load singularity 
    singularity exec ${params.image} actc -j \${PBS_NCPUS} ${subread} ${ccs} "${id}.subreads_to_ccs.bam"
    """

}

process deepconsensus
{

    input:
    tuple val (subread) , path (subreads_to_ccs) , path (ccs) , val  (id)

    output:
    tuple val (subread) , path ("${id}.output.fastq")

    script:
    """
    module load singularity 
    singularity exec ${params.image} deepconsensus run --subreads_to_ccs=${subreads_to_ccs}  --ccs_bam=${ccs} --checkpoint=/g/data/xl04/ka6418/temp/deepconsensus/testdata/model/checkpoint --cpus \${PBS_NCPUS} --output=${id}.output.fastq
    """

}


process concatFastq {

    input:
    tuple val (subread), path (fastqFiles)

    output:
    tuple val (subread) , path ("*.fastq")

    script:
    """
    subreadbase=\$(basename "${subread}" .bam)
    cat ${fastqFiles.join(' ')} > "\${subreadbase}.fastq"

    """
}

process cutadapt {

    publishDir "${params.output}", mode: 'copy'

    input: 
    tuple val (subread) , path (fastq)

    output:
    tuple val (subread), path ("*.trimmed.fastq.gz") , path ("*.md5sum") , path ("*.cutadapt.json")

    script:
    """
    module load cutadapt

    sample=\$(basename "${subread}" .bam)

    cutadapt --cores \${PBS_NCPUS} --anywhere file:\${PACBIOADAPTERS} \
    --error-rate 0.1 --overlap 25 --match-read-wildcards --revcomp --discard-trimmed \
    --json \${sample}.cutadapt.json \
    -o "\${sample}_deepconsensus.trimmed.fastq" \
    ${fastq}

    md5sum "\${sample}_deepconsensus.trimmed.fastq" > "\${sample}.md5sum"

    pigz "\${sample}_deepconsensus.trimmed.fastq"
 
    md5sum "\${sample}_deepconsensus.trimmed.fastq.gz" >> "\${sample}.md5sum"

    """

}

process qc {

    publishDir "${params.output}", mode: 'copy'

    input:
    tuple val (subread) , path (trimmed_fastq)

    output:
    tuple val (subread) , path ("*length_freq.csv") , path ("*quality_freq.csv") , path ("*stats.csv")

    script:

    """
    module load pythonlib

    sample=\$(basename "${subread}" .bam)
    python3 /g/data/te53/ncigflow/deepconsensus/bin/long_read_qc.py -input ${trimmed_fastq} -sample \${sample}
    """

}


workflow {

     
    subreadsChannel = Channel.from(params.subreads.split(':'))

    indexChannel = pbindex(subreadsChannel)
    
    indexChannel = indexChannel.map { tuple ->

                   key = groupKey(tuple[0],1)
                   
                   [key,tuple[1]]
                    }.groupTuple().flatMap { subreads_path, pbi_path ->
                                            pbi_path.collect {pbi ->
                                            [subreads_path, pbi]
                                            }
                        }

    idChannel = Channel.from(1..(params.chunks as Integer))
    indexChannel = indexChannel.combine(idChannel)

    
    ccs_bam = ccs(indexChannel)


    ccs_bam = ccs_bam.map {tuple ->

                key = groupKey(tuple[0],params.chunks)

               [key,tuple[1],tuple[2]]
    }.groupTuple().flatMap { it ->
            def (first, second, third) = it
            return [second, third].transpose().collect { [first, it[0], it[1]] }
        }

                            
    subreads_to_ccs = actc(ccs_bam)

    subreads_to_ccs = subreads_to_ccs.map {tuple ->

                key = groupKey(tuple[0],params.chunks)

               [key,tuple[1],tuple[2],tuple[3]]
    }.groupTuple().flatMap { it ->
            def (first, second, third, fourth) = it
            return [second, third, fourth].transpose().collect { [first, it[0], it[1], it[2]] }
        }


    fastq = deepconsensus(subreads_to_ccs)

    fastq = fastq.map {tuple ->

          key = groupKey(tuple[0],params.chunks)

          [key,tuple[1]]
    
    }.groupTuple()

    deepconsensus_fastq = concatFastq(fastq)

    deepconsensus_fastq = deepconsensus_fastq.map {tuple ->

          key = groupKey(tuple[0],1)

          [key,tuple[1]]
    
    }.groupTuple()

    cutadapt_fastq = cutadapt(deepconsensus_fastq)

    cutadapt_fastq = cutadapt_fastq.map {tuple ->

          key = groupKey(tuple[0],1)

          [key,tuple[1]]
    
    }.groupTuple()

    fastq_qc = qc(cutadapt_fastq)

}
