/*
 * Processes related to sequence quality assessment,
 *   quality control (e.g. trimming).
 */

process FASTQC {
  conda 'envs/fastqc.yml'
  cpus { read.size() } // number of read files to process
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

process MULTIQC {
  label "fastLocal"
  conda 'envs/multiqc.yml'
  publishDir "${params.publishDirPath}/multiQC/${params.multiQCLabel}"

  input:
    path(fastqc) // any number of fastqc files
  output:
    path("raw_multiqc_report.html")
  script:
    """
    multiqc -n raw_multiqc_report .
    """

}

process TRIMGALORE {
  conda 'envs/trim_galore.yml'
  label "fastLocal"

  input:
    tuple val(sample), path(forward_read), path(reverse_read)
  output:
    tuple val(sample), path("${ forward_read.simpleName }_val_1.fq.gz"), path("${ reverse_read.simpleName }_val_2.fq.gz"), path("${ forward_read }_trimming_report.txt"), path("${ reverse_read }_trimming_report.txt")
  script:
    /*
     * comments -> --ilumina # if adapters are not illumina, replace with adapters
     *   --paired  # only for PE studies, # if SE use only single read file
     */
    """
    trim_galore --gzip \
    --illumina \
    --phred33 \
    --paired $forward_read $reverse_read \
    --output_dir .
    """
}
