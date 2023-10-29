/*
 * This process is defined to clean up work files in a Nextflow pipeline.
 * It leverages an external bash script (clean_work_files.sh) to perform the cleanup.
 *
 * The process takes in a file as input and returns a value indicating the success of the cleanup.
 * Caching is set to 'lenient' which means the task will not fail if the cache is not available.
 * 
 * `baseDirPath` is expected to be provided as a parameter to the process, and it represents
 * the path where the clean_work_files.sh script resides.

 * Authors: Kosar Hooshmand
 */
 
import groovy.json.JsonBuilder

process clean_work_files {

  cache 'lenient'

  input:
  val(file)

  output:
  val(1), emit: IS_CLEAN

  script:
  """
    ${params.baseDirPath}/clean_work_files.sh "${file}"
  """
}

