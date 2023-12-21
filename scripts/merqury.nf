process meryl {

    executor = 'pbspro'
    queue = 'normal'
    project = 'xl04'
    time = '10h'
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
    outputfolder=${output}/"\$(basename "${file}")"_"${kmer}"
    /g/data/xl04/ka6418/bassiana/basequality/bassiana_illumina_notrim/meryl-1.4.1/bin/meryl count threads=\${PBS_NCPUS} k=${kmer} "${file}" output \${outputfolder}.meryl

    """

}

process meryl_unionsum {

    executor = 'pbspro'
    queue = 'normal'
    project = 'xl04'
    time = '1h'
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
    /g/data/xl04/ka6418/bassiana/basequality/bassiana_illumina_notrim/meryl-1.4.1/bin/meryl union-sum threads=\${PBS_NCPUS} output "${output}/${sample}_${tech}_${kmer}_combined.meryl" ${output}/*.meryl

    """

}

process merqury {

    executor = 'pbspro'
    queue = 'normal'
    project = 'xl04'
    time = '10h'
    clusterOptions = '-l ncpus=8,mem=16GB,storage=gdata/if89+gdata/xl04'

    input:
    val (fasta)
    val (dataset)
    val (output)
    val (sample)

    script:
    """
    # Activate conda environment
    source /g/data/xl04/ka6418/miniconda/etc/profile.d/conda.sh
    conda activate nanoplot_env

    # Run Merqury
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