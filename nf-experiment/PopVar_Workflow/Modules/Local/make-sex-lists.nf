process makesexlists {

    cache 'lenient'

    input:
    path gvcf_list 

    output:
    path("female.list"), optional: true, emit: female
    path("male.list"), optional: true, emit: male

    script:
    """
    ${params.baseDirPath}/make-sex-lists.sh $gvcf_list
    """
}
