process setup_directory {
        input:
        val topfolder 
        
        script:
        """
        mkdir -p ${topfolder}/assembly
        mkdir -p ${topfolder}/scaffolding
        mkdir -p ${topfolder}/scaffolding/yahs
        mkdir -p ${topfolder}/scaffolding/yahs_hicmap
        mkdir -p ${topfolder}/scaffolding/arima
        mkdir -p ${topfolder}/rawdata/shortread
        mkdir -p ${topfolder}/rawdata/shortread/trimmed
        mkdir -p ${topfolder}/rawdata/shortread/qc
        mkdir -p ${topfolder}/rawdata/shortread/plots
        mkdir -p ${topfolder}/rawdata/longread
        mkdir -p ${topfolder}/rawdata/longread/qc
        mkdir -p ${topfolder}/rawdata/longread/plots
        mkdir -p ${topfolder}/rawdata/kmers
        mkdir -p ${topfolder}/evaluation
        mkdir -p ${topfolder}/misc
        mkdir -p ${topfolder}/logs

        """
    }

workflow {

    setup_directory(params.topfolder)
}