nextflow.enable.dsl=2

include { DOWNLOAD_RAW_READS;
          DOWNLOAD_GENOME_ANNOTATIONS } from './modules/download.nf'
include { FASTQC as RAW_FASTQC } from './modules/quality.nf'
include { FASTQC as TRIM_FASTQC } from './modules/quality.nf'
include { MULTIQC as RAW_MULTIQC } from './modules/quality.nf' addParams(multiQCLabel: 'raw')
include { MULTIQC as TRIM_MULTIQC } from './modules/quality.nf' addParams(multiQCLabel: 'trimmed')
include { TRIMGALORE } from './modules/quality.nf'
include { BUILD_STAR } from './modules/genome.nf'

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

    GET_DATA.out | map{ it -> [ it[0], it[1][0], it[1][1] ] } | TRIMGALORE

    TRIMGALORE.out.reads | map{ it -> [ it[0], [ it[1], it[2] ] ]} | TRIM_FASTQC

    TRIM_FASTQC.out | map { it -> [ it[1], it[2] ] } \
                    | flatten \
                    | unique \
                    | collect \
                    | TRIM_MULTIQC \
                    | view

<<<<<<< HEAD
    DOWNLOAD_GENOME_ANNOTATIONS | view
=======
    DOWNLOAD_GENOME_ANNOTATIONS | BUILD_STAR | view
>>>>>>> 49175f75bc933d5d0ff07a14cd2ded933fc2c660

}
