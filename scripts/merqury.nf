process meryl {

    executor = 'pbspro'
    queue = 'normal'
    project = 'xl04'
    time = '48h'
    clusterOptions = '-l ncpus=48,mem=192GB,storage=gdata/if89+gdata/xl04'

    input:

    val (file)
    val (sample)
    val (tech)
    val (kmer)
    val (output)

    output:
    
    val "${output}/*.meryl"
    

    script:
    """
    module load merqury
    outputfolder=${output}/"\$(basename "${file}")"_"${kmer}"
    meryl count threads=\${PBS_NCPUS} k=${kmer} "${file}" output \${outputfolder}.meryl

    """

}

process meryl_unionsum {

    executor = 'pbspro'
    queue = 'normal'
    project = 'xl04'
    time = '48h'
    clusterOptions = '-l ncpus=48,mem=192GB,storage=gdata/if89+gdata/xl04'

    input:
    val (tocombine)
    val (sample)
    val (tech)
    val (kmer)
    val (output)
    

    output:
    val ("${output}/${sample}_${tech}_${kmer}_combined.meryl")

 

    script:
    """
    module load merqury
    meryl union-sum threads=\${PBS_NCPUS} output "${output}/${sample}_${tech}_${kmer}_combined.meryl" ${output}/*.meryl

    """

}

process merqury {

    executor = 'pbspro'
    queue = 'normal'
    project = 'xl04'
    time = '48h'
    clusterOptions = '-l ncpus=8,mem=16GB,storage=gdata/if89+gdata/xl04'

    input:
    val (fasta)
    val (dataset)
    val (output)
    val (sample)

    script:
    """

    module load merqury
    cd $output
    \${MERQURY}/merqury.sh ${dataset} ${fasta} ${sample}

    """

}


workflow {

    fileList = params.fileList.split(':')
    fileChannel = Channel.from(fileList)
    
    datasets = meryl(fileChannel,params.sample,params.tech,params.kmer,params.output)

    datasets.collect().map { it.join(' ') }.set { concatenatedDatasets }

    data = meryl_unionsum(concatenatedDatasets,params.sample,params.tech,params.kmer,params.output)
    merqury(params.fasta,data,params.output,params.sample)


}