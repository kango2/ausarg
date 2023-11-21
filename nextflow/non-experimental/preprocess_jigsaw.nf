process shortread_trimming {
    conda '/g/data/xl04/ka6418/miniconda/envs/genejigsaw'
    executor = 'pbspro'
    queue = 'express'
    project = 'xl04'
    time = '1h'
    clusterOptions = '-l ncpus=1,mem=2GB,storage=gdata/if89+gdata/xl04'

    input:
    tuple path (R1),path (R2),val (sample),val (flowcell),path (output),val (platform)

    output:
    tuple path ("$output/*_val_1.fq.gz"), path ("$output/*_val_2.fq.gz") , val (sample) , val (flowcell) , path (output) , val (platform) 


    script:

    """
    source /g/data/xl04/ka6418/miniconda/etc/profile.d/conda.sh
    conda activate genejigsaw
    trim_galore -o '${output}' --cores \${PBS_NCPUS} --paired '${R1}' '${R2}'

    """

}
