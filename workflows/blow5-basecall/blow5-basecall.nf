samplesch = Channel.from(params.files.split(':'))

models = ['HAC' , 'SUP']

process basecall {

    errorStrategy 'ignore' 
    
    publishDir "${params.outdir}" , mode: 'copy', pattern : "*fastq*"

    input:
    val (file)
    each mode

    output:
    path ("*${mode}*fastq.gz")

    script:

    """

    if [[ ${mode} == "HAC" ]]; then
        if [[ ${file} == *sqk-lsk109* || ${file} == *sqk-lsk110* || ${file} == *sqk-rad004* || ${file} == *sqk-ulk001* ]]; then
            MODEL="dna_r9.4.1_450bps_hac.cfg"  
        else
            MODEL="dna_r10.4.1_e8.2_400bps_hac.cfg"  
        fi

    elif [[ ${mode} == "SUP" ]]; then
        if [[ ${file} == *sqk-lsk109* || ${file} == *sqk-lsk110* || ${file} == *sqk-rad004* || ${file} == *sqk-ulk001* ]]; then
            MODEL="dna_r9.4.1_450bps_sup.cfg"  
        else
            MODEL="dna_r10.4.1_e8.2_400bps_5khz_sup.cfg"  
        fi
    fi

    export MODEL

    MODE=${mode}
    export MODE

    BASECALL_OUT=\${PWD}
    export BASECALL_OUT

    MERGED_SLOW5=${file}
    export MERGED_SLOW5

    bash /g/data/xl04/ka6418/github/ausarg/workflows/blow5-basecall/bin/basecall.sh

    """
}

process qc {

    publishDir "${params.outdir}" , mode: 'copy'

    input:
    val (fastq)

    output:
    path ("*.csv")

    script:

    """
    module load pythonlib
    python3 /g/data/xl04/ka6418/github/ausarg/workflows/blow5-basecall/bin/long_read_qc.py -input ${fastq} -output \${PWD} -sample \$(basename ${fastq} .fastq.gz)
  
    """

}


workflow {
    
    fastqCh = basecall(samplesch,models).flatten()
    qcCh = qc(fastqCh)
    


}