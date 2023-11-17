
include { fromQuery } from 'plugin/nf-sqldb'


process fastqmetrics {
    publishDir "/g/data/xl04/ka6418/github/ausarg/nextflow/pipeline", mode: 'copy', overwrite: false
    executor = 'pbspro'
    queue = 'normal'
    project = 'xl04'
    time = '1h'
    clusterOptions = '-l ncpus=8,mem=16GB,storage=gdata/if89+gdata/xl04'

    input:
    tuple path (fastq_file),val (sample),val (flowcell),path (output),val (platform)

    output:
    path ("$output/*_stats.csv")

    script:

    """
    module load biopython/1.79
    
    python3 /g/data/xl04/ka6418/github/ausarg/nextflow/long-read-qv/long-read-qv.py --i $fastq_file --o $output --s $sample --p $platform --f $flowcell

    """

}

workflow {
    

    fastq = channel
        .fromQuery('select title, flowcell, platform, filename from SRA', db: 'inputdb')
        .map { row ->
            // Assuming the row order is: title, flowcell, platform, filename
            def (title, flowcell, platform, filename) = row
            def output = "/g/data/xl04/ka6418/github/ausarg/nextflow/pipeline" // Define output path based on the title or any other unique identifier
            return [filename, title, flowcell, output, platform]
        }.view()

    fastqmetrics(fastq)
}




