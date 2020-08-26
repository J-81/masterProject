nextflow.enable.dsl=2

include { DOWNLOAD_RAW_READS } from './modules/download.nf'
include { FASTQC as RAW_FASTQC } from './modules/quality.nf'
include { FASTQC as TRIM_FASTQC } from './modules/quality.nf'
include { MULTIQC as RAW_MULTIQC } from './modules/quality.nf' addParams(multiQCLabel: 'raw')
include { MULTIQC as TRIM_MULTIQC } from './modules/quality.nf' addParams(multiQCLabel: 'trimmed')
include { TRIMGALORE } from './modules/quality.nf'

samples_ch = Channel.fromList( params.samples )
                    .take( params.limiter )

/*
 * Starting point, includes downloads data from GeneLab
 */

workflow GET_DATA {
  take: samples_ch
  main:
    samples_ch |  DOWNLOAD_RAW_READS
  emit:
    DOWNLOAD_RAW_READS.out.raw_reads
}

workflow {
	main:
    GET_DATA( samples_ch ) | RAW_FASTQC

    RAW_FASTQC.out | map { it -> [ it[1], it[2] ] } \
                   | flatten \
                   | unique \
                   | collect \
                   | RAW_MULTIQC \
                   | view

    GET_DATA.out | flatten | TRIMGALORE | TRIM_FASTQC

    TRIM_FASTQC.out | map { it -> [ it[1], it[2] ] } \
                    | flatten \
                    | unique \
                    | collect \
                    | TRIM_MULTIQC \
                    | view

}
