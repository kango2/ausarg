process longread_plots {

    executor = 'pbspro'
    queue = 'express'
    project = 'xl04'
    time = '1h'
    clusterOptions = '-l ncpus=1,mem=2GB,storage=gdata/if89+gdata/xl04'

    input:
    val (qc_folder)

    script:
    """
    module load Rlib

    Rscript /g/data/xl04/ka6418/github/ausarg/nextflow/automated_plotting.R -t ${qc_folder} -o /g/data/xl04/ka6418/nextflow_testing

    """


}

process kmer_plots {

    executor = 'pbspro'
    queue = 'express'
    project = 'xl04'
    time = '1h'
    clusterOptions = '-l ncpus=16,mem=192GB,storage=gdata/if89+gdata/xl04,jobfs=100GB'
    
    script:
    """
    /g/data/xl04/ka6418/github/ausarg/nextflow/kmer_nf.sh -i /g/data/xl04/ka6418/nextflow_testing/testhifi.fq.gz -s BASDU -o /g/data/xl04/ka6418/nextflow_testing -l 17 -t PacBio

    """


}

workflow {
    
    kmer_plots()

}