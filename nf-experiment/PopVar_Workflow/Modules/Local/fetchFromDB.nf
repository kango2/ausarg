import groovy.json.JsonBuilder

process fetchdataFromDB {

    cache 'lenient'

    input:
    val dummy

    output:
    path("output.tsv"), emit: dbOutput

    script:
    """
    python3 ${params.baseDirPath}/fetch_from_db.py > output.tsv
    """
}
