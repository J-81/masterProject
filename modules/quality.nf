/*
 * Processes related to sequence quality assessment,
 *   quality control (e.g. trimming).
 */

process FASTQC {
  conda 'envs/fastqc.yml'
  label 'cpuBound'

  input:
    tuple val(sample), path(raw_read)
  output:
    tuple val(sample), path("${raw_read.getSimpleName}.fastqc.html"), path("${raw_read.getSimpleName}.fastqc.zip")
  script:
    """
    fastqc -o . \
     -t $task.cpus \
      $raw_read
    """

}
