/*
 * Processes related to sequence quality assessment,
 *   quality control (e.g. trimming).
 */

process FASTQC {
  conda 'envs/fastqc.yml'
  label 'cpuBound'
  storeDir "${params.storeDirPath}/fastqc"

  input:
    tuple val(sample), path(read)
  output:
    tuple val(sample), path("*.fastqc.html"), path("*.fastqc.zip")
  script:
    """
    fastqc -o . \
     -t $task.cpus \
      $read
    """

}
