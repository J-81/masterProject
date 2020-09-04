nextflow.enable.dsl=2

include { DOWNLOAD_RAW_READS;
          DOWNLOAD_GENOME_ANNOTATIONS } from './modules/download.nf'
include { FASTQC as RAW_FASTQC } from './modules/quality.nf'
include { FASTQC as TRIM_FASTQC } from './modules/quality.nf'
include { MULTIQC as RAW_MULTIQC } from './modules/quality.nf' addParams(multiQCLabel: 'raw')
include { MULTIQC as TRIM_MULTIQC } from './modules/quality.nf' addParams(multiQCLabel: 'trimmed')
include { TRIMGALORE } from './modules/quality.nf'
include { BUILD_STAR;
          ALIGN_STAR;
          BUILD_RSEM;
          COUNT_ALIGNED } from './modules/genome.nf'
include { DGE_BY_DESEQ2 } from './modules/dge.nf'

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

    DOWNLOAD_GENOME_ANNOTATIONS | BUILD_STAR

    TRIMGALORE.out.reads | combine( BUILD_STAR.out ) | ALIGN_STAR

    ALIGN_STAR.out.transcriptomeMapping | view

    DOWNLOAD_GENOME_ANNOTATIONS.out | BUILD_RSEM

    ALIGN_STAR.out.transcriptomeMapping | combine( BUILD_RSEM.out ) | set { aligned_ch }
    aligned_ch | view
    aligned_ch | COUNT_ALIGNED

    COUNT_ALIGNED.out.countsPerGene | map { it[1] } | collect | set { rsem_ch }

    isa_ch = channel.fromPath("external/GLDS-104_metadata_GLDS-104-ISA.zip")
    organism_ch = channel.fromPath("external/organisms.csv")
    external_ch = isa_ch.mix(organism_ch)

    external_ch | combine( rsem_ch ) | set { for_dge_ch }
    for_dge_ch | view
    for_dge_ch | DGE_BY_DESEQ2

}
