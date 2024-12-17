include { fromQuery } from 'plugin/nf-sqldb'
include { fast5toblow5; basecall } from '/g/data/xl04/ka6418/github/ausarg/workflows/assembly/modules/processONT.nf'


workflow {

if ("${params.ONTmode}" == "fast5") {

    fast5Ch = Channel
    .fromQuery('select sample, flowcell, chemistry, filepath from ONT where format = "FAST5"', db: 'inputdb')
    .map { row ->
        def (sample, flowcell, chemistry, filepath) = row
        meta = [sample:sample]
        return [meta, filepath, chemistry, flowcell]
    }

    blow5Ch = fast5toblow5(fast5Ch)
    fastqCh = basecall(blow5Ch[0])
    fastqCh[0].view()

}
else if ("${params.ONTmode}" == "blow5") {

    blow5Ch = Channel
    .fromQuery('select sample, flowcell, chemistry, filepath from ONT where format = "BLOW5"', db: 'inputdb')
    .map { row ->
        def (sample, flowcell, chemistry, filepath) = row
        meta = [sample:sample]
        return [meta, filepath, chemistry, flowcell]
    }

    fastqCh = basecall(blow5Ch)
    fastqCh[0].view()

}
else if ("${params.ONTmode}" == "fastq") {

    fastqCh = Channel
    .fromQuery('select sample, flowcell, filepath from ONT where format = "FASTQ"', db: 'inputdb')
    .map { row ->
        def (sample, flowcell, filepath) = row
        meta = [sample:sample]
        return [meta, filepath, flowcell]
    }
    
    fastqCh.view()

}

}

