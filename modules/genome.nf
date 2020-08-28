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
STAR --runThreadN ${task.cpus} \
--runMode genomeGenerate \
--limitGenomeGenerateRAM 35000000000 \
--genomeSAindexNbases 14 \
--genomeDir STAR_REF \
--genomeFastaFiles ${ genomeFasta } \
--sjdbGTFfile ${ genomeGtf } \
--sjdbOverhang 149
    """

}
