//Add PacBio/Illumina/ONT in the names of output files and then merge just based on names


process align_ont {

    executor = 'pbspro'
    queue = 'normal'
    project = 'te53'
    time = '10h'
    clusterOptions = '-l ncpus=48,mem=192GB,storage=gdata/if89+gdata/xl04+gdata/te53'

    input:

    val (reads)
    val (ref)
    val (outputdir)

    output:
    val ("${outputdir}/${outputName}.sorted.ont.bam")


    script:
    def refBase = ref.replaceAll(/\.fasta$/, "")
    def readsBase = reads.replaceAll(/\.fastq\.gz$/, "")
    def outputName = "${refBase}_${readsBase}"


    """

    module load minimap2 samtools

    output=\$(basename "${ref}" .fasta)_\$(basename "${reads}" .fastq.gz)
    minimap2 -Y -K 2000M -t \${PBS_NCPUS} -ax map-ont ${ref} ${reads} | samtools sort -@ \${PBS_NCPUS} - --reference ${ref} -O BAM --write-index -o ${outputdir}/"\${output}.sorted.ont.bam"##idx##${outputdir}/"\${output}.sorted.ont.bam.bai"


    """

}

process align_pb {

    executor = 'pbspro'
    queue = 'normal'
    project = 'te53'
    time = '10h'
    clusterOptions = '-l ncpus=48,mem=192GB,storage=gdata/if89+gdata/xl04'

    input:

    val (reads)
    val (ref)
    val (outputdir)

    output:
    val ("${outputdir}/${outputName}.sorted.pb.bam")


    script:
    def refBase = ref.replaceAll(/\.fasta$/, "")
    def readsBase = reads.replaceAll(/\.fastq\.gz$/, "")
    def outputName = "${refBase}_${readsBase}"


    """

    module load minimap2 samtools

    output=\$(basename "${ref}" .fasta)_\$(basename "${reads}" .fastq.gz)
    minimap2 -Y -K 2000M -t \${PBS_NCPUS} -ax map-hifi ${ref} ${reads} | samtools sort -@ \${PBS_NCPUS} - --reference ${ref} -O BAM --write-index -o ${outputdir}/"\${output}.sorted.pb.bam"##idx##${outputdir}/"\${output}.sorted.pb.bam.bai"


    """

}

process bwa_index {

    executor = 'pbspro'
    queue = 'normal'
    project = 'te53'
    time = '10h'
    clusterOptions = '-l ncpus=48,mem=192GB,storage=gdata/if89+gdata/xl04'

    input:

    val (ref)

    output:

    val ("index.done")

    script:
    """
    module load bwa-mem2 
    bwa-mem2 index ${ref}
    touch index.done
    """

}

process align_illumina {

    executor = 'pbspro'
    queue = 'normal'
    project = 'te53'
    time = '10h'
    clusterOptions = '-l ncpus=48,mem=192GB,storage=gdata/if89+gdata/xl04'

    input:

    tuple val (r1), val (r2)
    val (ref)
    val (outputdir)
    val (index)

    output:
    val ("${outputdir}/${outputName}.sorted.illum.bam")


    script:
    def refBase = ref.replaceAll(/\.fasta$/, "")
    def r1Base = r1.replaceAll(/\.fastq\.gz$/, "")
    def r2Base = r2.replaceAll(/\.fastq\.gz$/, "")
    def outputName = "${refBase}_${r1Base}_${r2Base}"


    """

    module load bwa-mem2 samtools

    output=\$(basename "${ref}" .fasta)_\$(basename "${r1}" .fastq.gz)_\$(basename "${r2}" .fastq.gz)
    bwa-mem2 mem -t \${PBS_NCPUS} ${ref} ${r1} ${r2} > ${outputdir}/"\${output}.sam"
    samtools view -buS -@ \${PBS_NCPUS} ${outputdir}/"\${output}.sam" | samtools sort -@ \${PBS_NCPUS} - -o ${outputdir}/"\${output}.sorted.illum.bam" && samtools index ${outputdir}/"\${output}.sorted.illum.bam"


    """
}

process merge_illumina {

    executor = 'pbspro'
    queue = 'normal'
    project = 'te53'
    time = '10h'
    clusterOptions = '-l ncpus=48,mem=192GB,storage=gdata/if89+gdata/xl04'

    input:
    
    val (mergebam)
    val (sample)
    val (outputdir)

    output:

    val ("${outputdir}/${sample}.merged.illum.bam")

    script:
    """
    module load samtools
    
    samtools merge --write-index --threads \${PBS_NCPUS} -o "${outputdir}/${sample}.merged.illum.bam" ${outputdir}/*.sorted.illum.bam

    """

}

process merge_ont {

    executor = 'pbspro'
    queue = 'normal'
    project = 'te53'
    time = '10h'
    clusterOptions = '-l ncpus=48,mem=192GB,storage=gdata/if89+gdata/xl04+gdata/te53'

    input:
    
    val (mergebam)
    val (sample)
    val (outputdir)

    output:

    val ("${outputdir}/${sample}.merged.ont.bam")

    script:
    """
    module load samtools
    
    samtools merge --write-index --threads \${PBS_NCPUS} -o "${outputdir}/${sample}.merged.ont.bam" ${outputdir}/*.sorted.ont.bam

    """

}

process merge_pb {

    executor = 'pbspro'
    queue = 'normal'
    project = 'te53'
    time = '10h'
    clusterOptions = '-l ncpus=48,mem=192GB,storage=gdata/if89+gdata/xl04'

    input:
    
    val (mergebam)
    val (sample)
    val (outputdir)

    output:

    val ("${outputdir}/${sample}.merged.pb.bam")

    script:
    """
    module load samtools
    
    samtools merge --write-index --threads \${PBS_NCPUS} -o "${outputdir}/${sample}.merged.pb.bam" ${outputdir}/*.sorted.pb.bam

    """

}


process coverage {

    executor = 'pbspro'
    queue = 'normal'
    project = 'te53'
    time = '1h'
    clusterOptions = '-l ncpus=1,mem=4GB,storage=gdata/if89+gdata/xl04'

    input:

    val (bam)
    val (outputdir)

    output:

    val ("${outputdir}/*.depth.bed")

    script:

    """
    
    bam=${bam}
    export bam

    outdir=${outdir}
    export outdir

    bash /g/data/xl04/ka6418/github/ausarg/scripts/bam_to_bedcov.sh



    """


}






workflow {

   pbFiles = Channel.from(params.pbFiles.split(':')).view()
   ontFiles = Channel.from(params.ontFiles.split(':')).view()
   illuminaFiles = Channel.from(params.illuminaFiles.split(':'))
                .map { row ->
                    def (r1, r2) = row.split(';')
                    return [r1, r2]
                }.view()
                
   outdir = params.outdir
   sample = params.sample
   ref = params.ref

   index = bwa_index(ref)

   ontAligned = align_ont(ontFiles,ref,outdir)
   pbAligned = align_pb(pbFiles,ref,outdir)
   illumAligned = align_illumina(illuminaFiles,ref,outdir,index)

   ontAligned.collect().map { it.join(' ') }.set { ontConcat }
   pbAligned.collect().map { it.join(' ') }.set { pbConcat }
   illumAligned.collect().map { it.join(' ') }.set { illumConcat }

   ontBAM = merge_ont(ontConcat,sample,outdir)
   illumBAM = merge_illumina(illumConcat,sample,outdir)
   pbBAM = merge_pb(pbConcat,sample,outdir)

   sortedMergedBAM = ontBAM.mix(illumBAM, pbBAM)
   ontBAM.mix(illumBAM, pbBAM)
   //coverage(sortedMergedBAM,outdir)

}



