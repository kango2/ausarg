
process fastqeval {

    input:
    val (meta), val (fastq), val (flowcell)

    output:
    val (meta), path ("*stats.csv"), path ("*freq.csv")

    stub:

    """

    touch ${meta.sample}.${flowcell}.stats.csv
    touch ${meta.sample}.${flowcell}.freq.csv

    """


}