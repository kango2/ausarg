process hifiasm_assembly {
    executor = 'pbspro'
    queue = 'express'
    project = 'xl04'
    time = '1h'
    clusterOptions = '-l ncpus=1,mem=2GB,storage=gdata/if89+gdata/xl04'

    input: 
    tuple path(pacbio), path(ont), path(hic_R1), path(hic_R2), val(sample), path(output)

    output:
    path (*_HifiASM.fasta)

    script 
    """

    /g/data/xl04/ka6418/bassiana/hifiasm_bassiana/hifiasm/hifiasm -t \${PBS_NCPUS} -o "$output/$sample" --ul $ont --h1 $hic_R1 --h2 $hic_R2  $pacbio 

    awk '/^S/{print ">"\$2;print \$3}' "$output/$sample.asm.hic.p_ctg.gfa" > "$output/$sample_HifiASM.fasta"

    """

    
}
