/*
 * Dowload raw paired end reads from Genelab website, assumes a constant format on
 *   Genelab side
 */

process DOWNLOAD_RAW_READS {
  label 'networkBound'
  storeDir "${params.storeDirPath}/raw_reads"

  input:
    val(sample)
  output:
    tuple val(sample), path("${sample}_R?_raw.fastq.gz"), emit: raw_reads
  script:
    """
    wget --no-check-certificate --quiet \
    -O ${sample}_R1_raw.fastq.gz \
    ${params.GLDS_URL_PREFIX}${sample}_R1_raw.fastq.gz${params.GLDS_URL_SUFFIX}

    wget --no-check-certificate --quiet \
    -O ${sample}_R2_raw.fastq.gz \
    ${params.GLDS_URL_PREFIX}${sample}_R2_raw.fastq.gz${params.GLDS_URL_SUFFIX}
    """
}
