process telomeres {
    conda '/g/data/xl04/ka6418/miniconda/envs/genejigsaw'
    executor = 'pbspro'
    queue = 'express'
    project = 'xl04'
    time = '1h'
    clusterOptions = '-l ncpus=8,mem=16GB,storage=gdata/if89+gdata/xl04'

    input: 
    val (fasta)
    val (output)

    output:
    val ("${output}/*_Telomeres.csv")

    script:
    """
    export input=$fasta
    export output=$output
    export permatch=90
    export copies=100

    bash /g/data/xl04/ka6418/github/ausarg/scripts/find_telomeres.sh
    """

}

process centromeres {
    conda '/g/data/xl04/ka6418/miniconda/envs/genejigsaw'
    executor = 'pbspro'
    queue = 'express'
    project = 'xl04'
    time = '1h'
    clusterOptions = '-l ncpus=8,mem=16GB,storage=gdata/if89+gdata/xl04'

    input: 
    val (fasta)
    val (output)

    output:
    val ("${output}/Summary*.csv")

    script:
    """
    export inputfasta=$fasta
    export outputdir=$output

    bash /g/data/xl04/ka6418/github/ausarg/scripts/centromeres.sh
    """
}

process sequence_table {
    executor = 'pbspro'
    queue = 'express'
    project = 'xl04'
    time = '1h'
    clusterOptions = '-l ncpus=8,mem=16GB,storage=gdata/if89+gdata/xl04'

    input:
    val (fasta)
    val (output)

    output:
    val ("${output}/*_seqtable")

    script:
    """
    module load pythonlib 
    python3 /g/data/xl04/ka6418/github/ausarg/scripts/asm_to_sequencetable.py -fasta $fasta -outputdir $output

    """

}

process align_and_depth {

    executor = 'pbspro'
    queue = 'express'
    project = 'xl04'
    time = '1h'
    clusterOptions = '-l ncpus=8,mem=16GB,storage=gdata/if89+gdata/xl04'

    input:
    val (fasta)
    val (fastq)
    val (platform)
    val (output)

    output:
    val (${output}/"*sorted.bam.binned.depth.csv")

    script:
    """

    export platform=$platform
    export reference=$reference 
    export rawreads=$fastq
    export output=$output

    /g/data/xl04/ka6418/github/ausarg/scripts/align_and_depth.sh

    """

}

process gc {
    executor = 'pbspro'
    queue = 'express'
    project = 'xl04'
    time = '1h'
    clusterOptions = '-l ncpus=8,mem=16GB,storage=gdata/if89+gdata/xl04'

    input:
    val (fasta)
    val (output)

    output:
    val ("${output}/*_GC.csv")

    script:
    """
    export input=$fasta
    export output=$output

    bash /g/data/xl04/ka6418/github/ausarg/scripts/gc_content.sh
    """

}

workflow{
    //telomeres("/g/data/xl04/ka6418/nextflow_testing/testdata/BASDU_HifiASM.fasta","/g/data/xl04/ka6418/nextflow_testing/pipelinetest")
    //centromeres("/g/data/xl04/ka6418/nextflow_testing/testdata/BASDU_HifiASM.fasta","/g/data/xl04/ka6418/nextflow_testing/pipelinetest")
    sequence_table("/g/data/xl04/ka6418/nextflow_testing/testdata/BASDU_HifiASM.fasta","/g/data/xl04/ka6418/nextflow_testing/pipelinetest")
}