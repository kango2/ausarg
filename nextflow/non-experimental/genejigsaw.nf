include { fromQuery } from 'plugin/nf-sqldb'
include { longread_qc; shortread_qc} from '/g/data/xl04/ka6418/github/ausarg/nextflow/non-experimental/dataqc_jigsaw.nf'



workflow {
    
    longread_fastqs = channel
        .fromQuery('select title, flowcell, platform, filename from SRA where platform is "PACBIO_SMRT"', db: 'inputdb')
        .map { row ->
            def (title, flowcell, platform, filename) = row
            def output = "${params.topfolder}/longread_qc" 
            return [filename, title, flowcell, output, platform]
        }.view()

    longread_qc(longread_fastqs)


    illumina_fastqs = channel
    .fromQuery('select title, flowcell, platform, filename from SRA where platform is "ILLUMINA"', db: 'inputdb')
    .map { row ->
        def (title, flowcell, platform, filename) = row
        def output = "${params.topfolder}/shortread_qc"
        
        // Split the filename into two separate files
        def (file1, file2) = filename.split(':')
   
        // Return the split filenames along with other parameters
        return [file1, file2, title, flowcell, output, platform]
    }.view()

    shortread_qc(illumina_fastqs)


}


