params.image = "/g/data/xl04/ka6418/docker_images/deepconsensus_1.2.0-gpu.sif"
params.chunks = 10 

process pbindex {

    executor = 'pbspro'
    queue = 'normalsr'
    project = 'xl04'
    time = '24h'
    clusterOptions = '-l ncpus=104,mem=512GB,storage=gdata/if89+gdata/xl04'

    publishDir "${params.workdir}", mode:'copy'

    input:

    path (subreads)

    output:

    path ("${subreads}.pbi")

    script:

    """
    module load singularity 
    singularity exec ${params.image} pbindex ${subreads} -j \${PBS_NCPUS}

    """

}


process ccs {

    executor = 'pbspro'
    queue = 'normalsr'
    project = 'xl04'
    time = '24h'
    clusterOptions = '-l ncpus=104,mem=512GB,jobfs=400GB,storage=gdata/if89+gdata/xl04'

    publishDir "${params.output}", mode:'copy'

    input: 

    tuple path (pbi), path (subreads), val (id)

    output :

    path "${id}.ccs.bam"
    val (id)

    script:

    """
    module load singularity 
    singularity exec ${params.image} ccs --min-rq=0.88 -j \${PBS_NCPUS} --chunk="${id}"/"${params.chunks}" ${params.subreads} "${id}.ccs.bam"

    """


}


process actc 
{

    executor = 'pbspro'
    queue = 'normalsr'
    project = 'xl04'
    time = '24h'
    clusterOptions = '-l ncpus=104,mem=512GB,jobfs=400GB,storage=gdata/if89+gdata/xl04'

    publishDir "${params.output}", mode:'copy'

    input:
    path (ccs)
    val  (id)

    output:
    path "${id}.subreads_to_ccs.bam"
    path "${ccs}"
    val (id)

    script:
    """
    module load singularity 
    singularity exec ${params.image} actc -j \${PBS_NCPUS} ${params.subreads} ${ccs} "${id}.subreads_to_ccs.bam"
    """

}

process deepconsensus
{

    executor = 'pbspro'
    queue = 'gpuvolta'
    project = 'xl04'
    time = '24h'
    clusterOptions = '-l ncpus=12,mem=96GB,ngpus=1,jobfs=100GB,storage=gdata/if89+gdata/xl04'

    publishDir "${params.output}", mode:'copy'

    input:
    path (subreads_to_ccs)
    path (ccs)
    val  (id)

    output:
    path "${id}.output.fastq"

    script:
    """
    module load singularity 
    singularity exec ${params.image} deepconsensus run --subreads_to_ccs=${subreads_to_ccs}  --ccs_bam=${ccs} --checkpoint=/g/data/xl04/ka6418/temp/deepconsensus/testdata/model/checkpoint --cpus \${PBS_NCPUS} --output=${id}.output.fastq
    """

}


process concatFastq {

    executor = 'pbspro'
    queue = 'normal'
    project = 'xl04'
    time = '10h'
    clusterOptions = '-l ncpus=1,mem=4GB,storage=gdata/if89+gdata/xl04'

    publishDir "${params.output}", mode: 'copy'

    input:
    path fastqFiles
    val sample

    output:
    path "${sample}_deepconsensus.fastq"

    script:
    """
    cat ${fastqFiles.join(' ')} > "${sample}_deepconsensus.fastq"

    """
}

workflow {
    index = pbindex(params.subreads)
    idChannel = Channel.from(1..(params.chunks as Integer))
    combinedccs = index.combine(idChannel)
    ccs_bam = ccs(combinedccs.map { index, id -> [index, params.subreads, id] })
    subreads_to_ccs = actc(ccs_bam)
    fastq = deepconsensus(subreads_to_ccs)
    allfastq = fastq.collect()
    concatFastq(allfastq,params.sample)
}