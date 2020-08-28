/*
 * Processes related to genome and annotations
 */

process BUILD_STAR {
  conda 'envs/star.yml'
  label 'maxCPU'

  input:
    tuple path(genomeFasta), path(genomeGtf)
  output:
    path("STAR_REF/*")
  script:
    """
         # Number of available cores on server node
    STAR --runThreadN ${task.cpus} \
         --runMode genomeGenerate \
         # min needed for mouse ref
         --limitGenomeGenerateRAM 35000000000 \
         # this value is min(14,log2(GenomeLength)/2 - 1)
         # so for smaller genomes like bacteria this value will be smaller
         --genomeSAindexNbases 14 \
         --genomeDir STAR_REF \
         --genomeFastaFiles ${ genomeFasta } \
         --sjdbGTFfile ${ genomeGtf } \
         --sjdbOverhang ReadLength-1
    """

}
