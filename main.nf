nextflow.enable.dsl=2

include { DOWNLOAD_RAW_READS } from './modules/download.nf'
include { FASTQC; MULTIQC } from './modules/quality.nf'


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
    GET_DATA( samples_ch ) | FASTQC

    FASTQC.out | map { it -> [ it[1], it[2] ] } | flatten | collect | view

}
