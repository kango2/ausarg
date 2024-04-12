params.workdir = "/g/data/xl04/ka6418/temp/deepconsensus/testdata"
params.output = "/g/data/xl04/ka6418/temp/deepconsensus/nf"
params.subreads = "/g/data/xl04/ka6418/temp/deepconsensus/testdata/n1000.subreads.bam"

process pbindex {

    executor = 'pbspro'
    queue = 'express'
    project = 'xl04'
    time = '1h'
    clusterOptions = '-l ncpus=1,mem=4GB,storage=gdata/if89+gdata/xl04'

    publishDir "${params.workdir}", mode:'copy'

    input:

    path (subreads)

    output:

    path ("${subreads}.pbi")

    script:

    """
    module load singularity 
    singularity exec /g/data/xl04/ka6418/docker_images/deepconsensus_1.2.0-gpu.sif pbindex ${subreads} -j \${PBS_NCPUS}

    """

}


process ccs {

    executor = 'pbspro'
    queue = 'express'
    project = 'xl04'
    time = '1h'
    clusterOptions = '-l ncpus=48,mem=192GB,storage=gdata/if89+gdata/xl04'

    publishDir "${params.output}", mode:'copy'

    input: 

    tuple path (pbi), path (subreads), val (id)

    output :

    path "${id}.ccs.bam"
    val (id)

    script:

    """
    module load singularity 
    singularity exec /g/data/xl04/ka6418/docker_images/deepconsensus_1.2.0-gpu.sif ccs --min-rq=0.88 -j \${PBS_NCPUS} --chunk="${id}"/"${params.chunks}" n1000.subreads.bam "${id}.ccs.bam"

    """


}


process actc 
{
    executor = 'pbspro'
    queue = 'express'
    project = 'xl04'
    time = '1h'
    clusterOptions = '-l ncpus=48,mem=192GB,storage=gdata/if89+gdata/xl04'

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
    singularity exec /g/data/xl04/ka6418/docker_images/deepconsensus_1.2.0-gpu.sif actc -j \${PBS_NCPUS} ${params.subreads} ${ccs} "${id}.subreads_to_ccs.bam"
    """

}

process deepconsensus
{
    executor = 'pbspro'
    queue = 'express'
    project = 'xl04'
    time = '1h'
    clusterOptions = '-l ncpus=48,mem=192GB,storage=gdata/if89+gdata/xl04'

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
    singularity exec /g/data/xl04/ka6418/docker_images/deepconsensus_1.2.0-gpu.sif deepconsensus run --subreads_to_ccs=${subreads_to_ccs}  --ccs_bam=${ccs} --checkpoint=/g/data/xl04/ka6418/temp/deepconsensus/testdata/model/checkpoint --output=${id}.output.fastq
    """

}


process concatFastq {

    executor = 'pbspro'
    queue = 'express'
    project = 'xl04'
    time = '1h'
    clusterOptions = '-l ncpus=48,mem=192GB,storage=gdata/if89+gdata/xl04'

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

    params.chunks = 10 
    idChannel = Channel.from(1..(params.chunks as Integer))
    combinedccs = index.combine(idChannel)

    ccs_bam = ccs(combinedccs.map { index, id -> [index, params.subreads, id] })
    subreads_to_ccs = actc(ccs_bam)

    fastq = deepconsensus(subreads_to_ccs)

    allfastq = fastq.collect()
    concatFastq(allfastq,"BASDU")

 

}