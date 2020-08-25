/*
 * Processes related to sequence quality assessment,
 *   quality control (e.g. trimming).
 */

process FASTQC {
  conda 'envs/fastqc.yml'
  cpus { read.length } // number of read files to process
  storeDir "${params.storeDirPath}/fastqc"

  input:
    tuple val(sample), path(read)
  output:
    tuple val(sample), path("*_fastqc.html"), path("*_fastqc.zip")
  script:
    """
    fastqc -o . \
     -t $task.cpus \
      $read
    """

}
