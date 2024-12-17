

process fast5toblow5 {

    input:
    tuple val (meta), val (fast5), val (chemistry), val (flowcell)

    output:
    tuple val (meta), path ("*.blow5"), val (chemistry), val (flowcell)
    tuple val (meta), path ("*.blow5.idx"), path ("*.blow5.md5")

    stub:

    """
    
    touch ${meta.sample}.${flowcell}.blow5
    touch ${meta.sample}.${flowcell}.blow5.idx
    touch ${meta.sample}.${flowcell}.blow5.md5

    """

}


process basecall {

    input:
    tuple val (meta), val (blow5) , val (chemistry), val (flowcell)

    output:
    tuple val (meta), path ("*pass.fastq.gz"), val (flowcell)
    tuple val (meta), path ("*.gzi"), path ("*fail*"), path ("*md5*")

    stub:

    """

    touch ${meta.sample}.${flowcell}.pass.fastq.gz
    touch ${meta.sample}.${flowcell}.pass.fastq.gz.gzi
    touch ${meta.sample}.${flowcell}.fail.fastq.gz
    touch ${meta.sample}.${flowcell}.fail.fastq.gz.gzi
    touch ${meta.sample}.${flowcell}.fastq.md5

    """

}